#include "restore.h"

#include "cmp.h"
#include "error.h"
#include "fail.h"
#include "log.h"
#include "mem.h"
#include "util.h"
#include "working_set.h"

#include "preprocessor/preprocessing.h"

struct SleqpRestoration
{
  int refcount;

  SleqpPreprocessingState* preprocessing_state;
  SleqpProblem* original_problem;
  SleqpProblem* transformed_problem;

  SleqpParams* params;

  SLEQP_ACTIVE_STATE* working_var_states;
  SLEQP_ACTIVE_STATE* working_cons_states;

  double* var_dual;
  double* cons_dual;

  double* cache;
  SleqpSparseVec* stationarity_residuals;
  double* dense_stationarity_residuals;

  // Maps from transformed to original
  int* linear_cons_map;
  int* var_map;
};

static SLEQP_RETCODE
create_maps(SleqpRestoration* restoration)
{
  SleqpPreprocessingState* state = restoration->preprocessing_state;

  SleqpProblem* problem = sleqp_preprocessing_state_get_problem(state);

  const int num_variables = sleqp_problem_num_vars(problem);
  const int num_linear    = sleqp_problem_num_lin_cons(problem);

  {
    SleqpVariableState* var_states
      = sleqp_preprocessing_state_variable_states(state);

    int offset = 0;

    for (int j = 0; j < num_variables; ++j)
    {
      if (var_states[j].state == SLEQP_VAR_UNCHANGED)
      {
        restoration->var_map[j - offset] = j;
      }
      else
      {
        ++offset;
      }
    }
  }

  {
    SleqpConstraintState* linear_cons_states
      = sleqp_preprocessing_state_linear_constraint_states(state);

    int offset = 0;

    for (int i = 0; i < num_linear; ++i)
    {
      if (linear_cons_states[i].state == SLEQP_CONS_UNCHANGED)
      {
        restoration->linear_cons_map[i - offset] = i;
      }
      else
      {
        ++offset;
      }
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_restoration_create(SleqpRestoration** star,
                         SleqpPreprocessingState* preprocessing_state,
                         SleqpProblem* transformed_problem,
                         SleqpParams* params)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpRestoration* restoration = *star;

  *restoration          = (SleqpRestoration){0};
  restoration->refcount = 1;

  restoration->preprocessing_state = preprocessing_state;
  SLEQP_CALL(
    sleqp_preprocessing_state_capture(restoration->preprocessing_state));

  SleqpProblem* problem
    = sleqp_preprocessing_state_get_problem(preprocessing_state);

  restoration->original_problem = problem;
  SLEQP_CALL(sleqp_problem_capture(restoration->original_problem));

  restoration->transformed_problem = problem;
  SLEQP_CALL(sleqp_problem_capture(restoration->transformed_problem));

  restoration->params = params;
  SLEQP_CALL(sleqp_params_capture(restoration->params));

  const int num_variables   = sleqp_problem_num_vars(problem);
  const int num_constraints = sleqp_problem_num_cons(problem);

  SLEQP_CALL(
    sleqp_alloc_array(&restoration->working_var_states, num_variables));
  SLEQP_CALL(
    sleqp_alloc_array(&restoration->working_cons_states, num_constraints));

  SLEQP_CALL(sleqp_alloc_array(&restoration->var_dual, num_variables));
  SLEQP_CALL(sleqp_alloc_array(&restoration->cons_dual, num_constraints));

  SLEQP_CALL(sleqp_alloc_array(&restoration->cache, num_variables));

  SLEQP_CALL(
    sleqp_sparse_vector_create_empty(&restoration->stationarity_residuals,
                                     num_variables));

  SLEQP_CALL(sleqp_alloc_array(&restoration->dense_stationarity_residuals,
                               num_variables));

  {
    const int num_variables = sleqp_problem_num_vars(transformed_problem);
    const int num_linear    = sleqp_problem_num_lin_cons(transformed_problem);

    SLEQP_CALL(sleqp_alloc_array(&restoration->linear_cons_map, num_linear));
    SLEQP_CALL(sleqp_alloc_array(&restoration->var_map, num_variables));

    SLEQP_CALL(create_maps(restoration));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
restore_primal(SleqpRestoration* restoration,
               const SleqpSparseVec* source,
               SleqpSparseVec* target)
{
  int num_fixed_vars;
  int* fixed_var_indices;
  double* fixed_var_values;

  SLEQP_CALL(
    sleqp_preprocessing_state_fixed_variables(restoration->preprocessing_state,
                                              &num_fixed_vars,
                                              &fixed_var_indices,
                                              &fixed_var_values));

  SLEQP_CALL(sleqp_preprocessing_merge_entries(source,
                                               target,
                                               num_fixed_vars,
                                               fixed_var_indices,
                                               fixed_var_values));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
store_working_set(const SleqpRestoration* restoration,
                  SleqpWorkingSet* working_set)
{
  SleqpProblem* problem = restoration->original_problem;

  const int num_variables   = sleqp_problem_num_vars(problem);
  const int num_constraints = sleqp_problem_num_cons(problem);

  SLEQP_CALL(sleqp_working_set_reset(working_set));

  for (int j = 0; j < num_variables; ++j)
  {
    SLEQP_ACTIVE_STATE state = restoration->working_var_states[j];

    if (state != SLEQP_INACTIVE)
    {
      SLEQP_CALL(sleqp_working_set_add_var(working_set, j, state));
    }
  }

  for (int i = 0; i < num_constraints; ++i)
  {
    SLEQP_ACTIVE_STATE state = restoration->working_cons_states[i];

    if (state != SLEQP_INACTIVE)
    {
      SLEQP_CALL(sleqp_working_set_add_cons(working_set, i, state));
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
store_duals(const SleqpRestoration* restoration, SleqpIterate* original)
{
  SleqpProblem* problem = restoration->original_problem;

  const int num_variables   = sleqp_problem_num_vars(problem);
  const int num_constraints = sleqp_problem_num_cons(problem);

  const double zero_eps
    = sleqp_params_value(restoration->params, SLEQP_PARAM_ZERO_EPS);

  SLEQP_CALL(sleqp_sparse_vector_from_raw(sleqp_iterate_vars_dual(original),
                                          restoration->var_dual,
                                          num_variables,
                                          zero_eps));

  SLEQP_CALL(sleqp_sparse_vector_from_raw(sleqp_iterate_cons_dual(original),
                                          restoration->cons_dual,
                                          num_constraints,
                                          zero_eps));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
reset_working_set_states(SleqpRestoration* restoration)
{
  SleqpProblem* problem = restoration->original_problem;

  const int num_variables   = sleqp_problem_num_vars(problem);
  const int num_constraints = sleqp_problem_num_cons(problem);

  for (int i = 0; i < num_constraints; ++i)
  {
    restoration->working_cons_states[i] = SLEQP_INACTIVE;
  }

  for (int j = 0; j < num_variables; ++j)
  {
    restoration->working_var_states[j] = SLEQP_INACTIVE;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
prepare_working_set(SleqpRestoration* restoration,
                    const SleqpWorkingSet* transformed,
                    SleqpWorkingSet* original)
{
  SleqpPreprocessingState* preprocessing_state
    = restoration->preprocessing_state;
  SleqpProblem* problem = restoration->original_problem;

  const int num_variables = sleqp_problem_num_vars(problem);

  const int num_general = sleqp_problem_num_gen_cons(problem);
  const int num_linear  = sleqp_problem_num_lin_cons(problem);

  SleqpConvertedBound* converted_bounds;
  int num_converted_bounds;

  SLEQP_CALL(sleqp_preprocessing_state_converted_bounds(preprocessing_state,
                                                        &converted_bounds,
                                                        &num_converted_bounds));

  SLEQP_CALL(reset_working_set_states(restoration));

  SleqpVariableState* var_states
    = sleqp_preprocessing_state_variable_states(preprocessing_state);

  SleqpConstraintState* linear_cons_states
    = sleqp_preprocessing_state_linear_constraint_states(preprocessing_state);

  {

    int offset = 0;

    for (int j = 0; j < num_variables; ++j)
    {
      if (var_states[j].state == SLEQP_VAR_UNCHANGED)
      {
        restoration->working_var_states[j]
          = sleqp_working_set_var_state(transformed, j - offset);
      }
      else
      {
        ++offset;
      }
    }
  }

  {
    for (int i = 0; i < num_general; ++i)
    {
      restoration->working_cons_states[i]
        = sleqp_working_set_cons_state(transformed, i);
    }

    int offset = 0;

    for (int i = 0; i < num_linear; ++i)
    {
      const int i_general = i + num_general;

      if (linear_cons_states[i].state == SLEQP_CONS_UNCHANGED)
      {
        restoration->working_cons_states[i_general]
          = sleqp_working_set_cons_state(transformed, i_general - offset);
      }
      else
      {
        ++offset;
      }
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
reset_duals(SleqpRestoration* restoration)
{
  SleqpProblem* problem = restoration->original_problem;

  const int num_variables   = sleqp_problem_num_vars(problem);
  const int num_constraints = sleqp_problem_num_cons(problem);

  for (int j = 0; j < num_variables; ++j)
  {
    restoration->var_dual[j] = SLEQP_INACTIVE;
  }

  for (int i = 0; i < num_constraints; ++i)
  {
    restoration->cons_dual[i] = SLEQP_INACTIVE;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
correct_variable_stationarity_residuum(SleqpRestoration* restoration,
                                       int variable)
{
  double* residuals = restoration->dense_stationarity_residuals;

  assert(restoration->var_dual[variable] == 0.);
  restoration->var_dual[variable] = -residuals[variable];
  residuals[variable]             = 0.;

  return SLEQP_OKAY;
}

static SLEQP_ACTIVE_STATE
desired_var_state_in_forcing_constraint(
  SleqpForcingConstraint* forcing_constraint,
  int k)
{
  const SleqpBoundState bound_state = forcing_constraint->state;

  const double factor = forcing_constraint->factors[k];

  const bool cons_at_upper = (bound_state == SLEQP_UPPER_BOUND);

  const bool prod_at_upper = !(cons_at_upper);

  const bool pos_factor = (factor > 0);

  const bool var_at_upper = pos_factor ? (prod_at_upper) : (!prod_at_upper);

  return var_at_upper ? SLEQP_ACTIVE_UPPER : SLEQP_ACTIVE_LOWER;
}

static SLEQP_RETCODE
correct_forcing_constraint(SleqpRestoration* restoration,
                           SleqpForcingConstraint* forcing_constraint)
{
  const int num_variables = forcing_constraint->num_variables;

  double* residuals = restoration->dense_stationarity_residuals;

  const double eps = sleqp_params_value(restoration->params, SLEQP_PARAM_EPS);

  SLEQP_NUM_ASSERT_PARAM(eps);

  bool has_inverted_residual_sign = false;
  int max_index;
  double max_dual = 0.;

  for (int k = 0; k < num_variables; ++k)
  {
    const int j = forcing_constraint->variables[k];

    SLEQP_ACTIVE_STATE desired_var_state
      = desired_var_state_in_forcing_constraint(forcing_constraint, k);

    const bool var_at_upper = (desired_var_state == SLEQP_ACTIVE_UPPER);

    const double target_dual = -residuals[j];
    const bool nonneg_dual   = (target_dual >= 0.);
    const bool correct_sign  = (var_at_upper == nonneg_dual);

    if (!correct_sign)
    {
      const double factor = forcing_constraint->factors[k];

      has_inverted_residual_sign = true;

      double current_dual = -residuals[j] / factor;

      if (SLEQP_ABS(current_dual) > SLEQP_ABS((max_dual)))
      {
        max_index = k;
        max_dual  = current_dual;
      }
    }
  }

  if (has_inverted_residual_sign)
  {
    const SleqpBoundState bound_state = forcing_constraint->state;

    SLEQP_ACTIVE_STATE cons_state;

    if (bound_state == SLEQP_UPPER_BOUND)
    {
      assert(max_dual >= 0.);
      cons_state = SLEQP_ACTIVE_UPPER;
    }
    else
    {
      assert(max_dual <= 0.);
      cons_state = SLEQP_ACTIVE_LOWER;
    }

    restoration->working_cons_states[forcing_constraint->constraint]
      = cons_state;

    const double cons_dual = max_dual;

    assert(restoration->cons_dual[forcing_constraint->constraint] == 0.);

    restoration->cons_dual[forcing_constraint->constraint] = cons_dual;

    for (int k = 0; k < num_variables; ++k)
    {
      const int j = forcing_constraint->variables[k];

      if (k == max_index)
      {
        residuals[j] = 0.;
        continue;
      }

      SLEQP_ACTIVE_STATE var_state
        = desired_var_state_in_forcing_constraint(forcing_constraint, k);

      restoration->working_var_states[j] = var_state;

      restoration->var_dual[j]
        = -(residuals[j] + forcing_constraint->factors[k] * cons_dual);

      if (var_state == SLEQP_ACTIVE_LOWER)
      {
        sleqp_assert_is_leq(restoration->var_dual[j], 0., eps);
      }
      else if (var_state == SLEQP_ACTIVE_UPPER)
      {
        sleqp_assert_is_geq(restoration->var_dual[j], 0., eps);
      }

      residuals[j] = 0.;
    }
  }
  else
  {
    for (int k = 0; k < num_variables; ++k)
    {
      const int j = forcing_constraint->variables[k];

      SLEQP_ACTIVE_STATE var_state
        = desired_var_state_in_forcing_constraint(forcing_constraint, k);

      restoration->working_var_states[j] = var_state;

      SLEQP_CALL(correct_variable_stationarity_residuum(restoration, j));
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
correct_converted_bound(SleqpRestoration* restoration,
                        SleqpConvertedBound* converted_bound)
{
  SleqpProblem* problem = restoration->original_problem;

  const int num_general = sleqp_problem_num_gen_cons(problem);

  const int i_general = converted_bound->constraint + num_general;

  const int j = converted_bound->variable;

  SLEQP_ACTIVE_STATE var_state = restoration->working_var_states[j];

  const bool bound_flip = converted_bound->factor < 0.;

  const double var_dual  = restoration->var_dual[j];
  const double cons_dual = var_dual / converted_bound->factor;

  assert(restoration->working_cons_states[i_general] == SLEQP_INACTIVE);
  assert(restoration->cons_dual[i_general] == 0.);

  if (var_state == SLEQP_ACTIVE_BOTH)
  {
    var_state = (var_dual >= 0) ? SLEQP_ACTIVE_UPPER : SLEQP_ACTIVE_LOWER;
  }

  if (var_state == SLEQP_INACTIVE)
  {
    assert(restoration->var_dual[j] == 0.);

    return SLEQP_OKAY;
  }
  else if (var_state == SLEQP_ACTIVE_LOWER)
  {
    if (converted_bound->state & SLEQP_LOWER_BOUND)
    {
      if (bound_flip)
      {
        restoration->working_cons_states[i_general] = SLEQP_ACTIVE_UPPER;
      }
      else
      {
        restoration->working_cons_states[i_general] = SLEQP_ACTIVE_LOWER;
      }

      restoration->cons_dual[i_general] = cons_dual;

      restoration->var_dual[j]           = 0.;
      restoration->working_var_states[j] = SLEQP_INACTIVE;
    }
  }
  else
  {
    assert(var_state == SLEQP_ACTIVE_UPPER);

    if (converted_bound->state & SLEQP_UPPER_BOUND)
    {
      if (bound_flip)
      {
        restoration->working_cons_states[i_general] = SLEQP_ACTIVE_LOWER;
      }
      else
      {
        restoration->working_cons_states[i_general] = SLEQP_ACTIVE_UPPER;
      }

      restoration->cons_dual[i_general] = cons_dual;

      restoration->var_dual[j]           = 0.;
      restoration->working_var_states[j] = SLEQP_INACTIVE;
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
correct_converted_bounds(SleqpRestoration* restoration)
{
  SleqpPreprocessingState* preprocessing_state
    = restoration->preprocessing_state;

  SleqpConvertedBound* converted_bounds;
  int num_converted_bounds;

  SLEQP_CALL(sleqp_preprocessing_state_converted_bounds(preprocessing_state,
                                                        &converted_bounds,
                                                        &num_converted_bounds));

  for (int k = 0; k < num_converted_bounds; ++k)
  {
    SLEQP_CALL(correct_converted_bound(restoration, converted_bounds + k));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
correct_forcing_constraints(SleqpRestoration* restoration)
{
  SleqpPreprocessingState* preprocessing_state
    = restoration->preprocessing_state;

  SleqpForcingConstraint* forcing_constraints;
  int num_forcing_constraints;

  SLEQP_CALL(
    sleqp_preprocessing_state_forcing_constraints(preprocessing_state,
                                                  &forcing_constraints,
                                                  &num_forcing_constraints));

  for (int k = 0; k < num_forcing_constraints; ++k)
  {
    SLEQP_CALL(
      correct_forcing_constraint(restoration, forcing_constraints + k));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
correct_fixed_variable(SleqpRestoration* restoration, int j)
{
  restoration->working_var_states[j] = SLEQP_ACTIVE_BOTH;

  SLEQP_CALL(correct_variable_stationarity_residuum(restoration, j));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
correct_fixed_variables(SleqpRestoration* restoration)
{
  int num_fixed_vars;
  int* fixed_var_indices;
  double* fixed_var_values;

  SleqpPreprocessingState* preprocessing_state
    = restoration->preprocessing_state;

  SLEQP_CALL(sleqp_preprocessing_state_fixed_variables(preprocessing_state,
                                                       &num_fixed_vars,
                                                       &fixed_var_indices,
                                                       &fixed_var_values));

  SleqpVariableState* var_states
    = sleqp_preprocessing_state_variable_states(preprocessing_state);

  for (int k = 0; k < num_fixed_vars; ++k)
  {
    const int j = fixed_var_indices[k];

    if (var_states[j].state == SLEQP_VAR_BOUND_FIXED)
    {
      SLEQP_CALL(correct_fixed_variable(restoration, j));
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
prepare_duals(SleqpRestoration* restoration, const SleqpIterate* transformed)
{
  SleqpProblem* problem = restoration->original_problem;

  const int num_general = sleqp_problem_num_gen_cons(problem);

  SLEQP_CALL(reset_duals(restoration));

  {
    SleqpSparseVec* var_dual = sleqp_iterate_vars_dual(transformed);

    for (int k = 0; k < var_dual->nnz; ++k)
    {
      int j    = var_dual->indices[k];
      double v = var_dual->data[k];

      restoration->var_dual[restoration->var_map[j]] = v;
    }
  }

  {
    SleqpSparseVec* cons_dual = sleqp_iterate_cons_dual(transformed);

    for (int k = 0; k < cons_dual->nnz; ++k)
    {
      int i    = cons_dual->indices[k];
      double v = cons_dual->data[k];

      if (i < num_general)
      {
        restoration->cons_dual[i] = v;
        continue;
      }

      i -= num_general;

      restoration->cons_dual[num_general + restoration->linear_cons_map[i]] = v;
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_stationarity_residuals(SleqpRestoration* restoration,
                               SleqpIterate* original)
{
  SleqpProblem* problem = restoration->original_problem;

  const double zero_eps
    = sleqp_params_value(restoration->params, SLEQP_PARAM_ZERO_EPS);

  SLEQP_CALL(store_duals(restoration, original));

  SleqpSparseVec* residuals = restoration->stationarity_residuals;

  SLEQP_CALL(sleqp_iterate_stationarity_residuals(problem,
                                                  original,
                                                  restoration->cache,
                                                  residuals,
                                                  zero_eps));

  double* dense_residuals = restoration->dense_stationarity_residuals;

  SLEQP_CALL(sleqp_sparse_vector_to_raw(residuals, dense_residuals));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_restoration_restore_iterate(SleqpRestoration* restoration,
                                  const SleqpIterate* transformed,
                                  SleqpIterate* original)
{
  SleqpProblem* problem = restoration->original_problem;

  SLEQP_CALL(restore_primal(restoration,
                            sleqp_iterate_primal(transformed),
                            sleqp_iterate_primal(original)));

  bool reject = false;

  SLEQP_CALL(sleqp_set_and_evaluate(problem,
                                    original,
                                    SLEQP_VALUE_REASON_NONE,
                                    &reject));

  if (reject)
  {
    sleqp_raise(SLEQP_FUNC_EVAL_ERROR, "Function rejected restored solution");
  }

  SLEQP_CALL(prepare_working_set(restoration,
                                 sleqp_iterate_working_set(transformed),
                                 sleqp_iterate_working_set(original)));

  SLEQP_CALL(prepare_duals(restoration, transformed));

  SLEQP_CALL(compute_stationarity_residuals(restoration, original));

  SLEQP_CALL(correct_fixed_variables(restoration));

  SLEQP_CALL(correct_converted_bounds(restoration));

  SLEQP_CALL(correct_forcing_constraints(restoration));

  SLEQP_CALL(
    store_working_set(restoration, sleqp_iterate_working_set(original)));

  SLEQP_CALL(store_duals(restoration, original));

#ifndef DEBUG

  assert(sleqp_working_set_valid(sleqp_iterate_working_set(original)));

  // This is only true if the working set of the transformed iterate is
  // supporting the respective duals as well.
  /*
  bool supports_cons_dual;
  SLEQP_CALL(sleqp_working_set_supports_cons_dual(sleqp_iterate_get_working_set(original),
                                                  sleqp_iterate_get_cons_dual(original),
                                                  &supports_cons_dual));

  assert(supports_cons_dual);

  bool supports_vars_dual;
  SLEQP_CALL(sleqp_working_set_supports_vars_dual(sleqp_iterate_get_working_set(original),
                                                  sleqp_iterate_get_vars_dual(original),
                                                  &supports_vars_dual));

  assert(supports_vars_dual);
  */

#endif

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_restoration_capture(SleqpRestoration* restoration)
{
  ++restoration->refcount;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
restoration_free(SleqpRestoration** star)
{
  SleqpRestoration* restoration = *star;

  sleqp_free(&restoration->var_map);
  sleqp_free(&restoration->linear_cons_map);

  sleqp_free(&restoration->dense_stationarity_residuals);

  SLEQP_CALL(sleqp_sparse_vector_free(&restoration->stationarity_residuals));

  sleqp_free(&restoration->cache);

  sleqp_free(&restoration->cons_dual);
  sleqp_free(&restoration->var_dual);

  sleqp_free(&restoration->working_cons_states);
  sleqp_free(&restoration->working_var_states);

  SLEQP_CALL(sleqp_params_release(&restoration->params));

  SLEQP_CALL(sleqp_problem_release(&restoration->transformed_problem));
  SLEQP_CALL(sleqp_problem_release(&restoration->original_problem));

  SLEQP_CALL(
    sleqp_preprocessing_state_release(&restoration->preprocessing_state));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_restoration_release(SleqpRestoration** star)
{
  SleqpRestoration* restoration = *star;

  if (!restoration)
  {
    return SLEQP_OKAY;
  }

  if (--restoration->refcount == 0)
  {
    SLEQP_CALL(restoration_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
