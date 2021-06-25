#include "sleqp_restore.h"

#include "sleqp_mem.h"
#include "sleqp_util.h"

#include "preprocessor/sleqp_preprocessing.h"

struct SleqpRestoration
{
  int refcount;

  SleqpPreprocessingState* preprocessing_state;
  SleqpProblem* original_problem;

  SleqpParams* params;

  SLEQP_ACTIVE_STATE* working_var_states;
  SLEQP_ACTIVE_STATE* working_cons_states;

  double* var_dual;
  double* cons_dual;

  double* cache;
  SleqpSparseVec* stationarity_residuals;
};


SLEQP_RETCODE sleqp_restoration_create(SleqpRestoration** star,
                                       SleqpPreprocessingState* preprocessing_state,
                                       SleqpParams* params)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpRestoration* restoration = *star;

  *restoration = (SleqpRestoration) {0};
  restoration->refcount = 1;

  restoration->preprocessing_state = preprocessing_state;
  SLEQP_CALL(sleqp_preprocessing_state_capture(restoration->preprocessing_state));

  SleqpProblem* problem = sleqp_preprocessing_state_get_problem(preprocessing_state);

  restoration->original_problem = problem;
  SLEQP_CALL(sleqp_problem_capture(restoration->original_problem));

  restoration->params = params;
  SLEQP_CALL(sleqp_params_capture(restoration->params));

  const int num_variables = sleqp_problem_num_variables(problem);
  const int num_constraints = sleqp_problem_num_constraints(problem);

  SLEQP_CALL(sleqp_alloc_array(&restoration->working_var_states, num_variables));
  SLEQP_CALL(sleqp_alloc_array(&restoration->working_cons_states, num_constraints));

  SLEQP_CALL(sleqp_alloc_array(&restoration->var_dual, num_variables));
  SLEQP_CALL(sleqp_alloc_array(&restoration->cons_dual, num_constraints));

  SLEQP_CALL(sleqp_alloc_array(&restoration->cache, num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&restoration->stationarity_residuals, num_variables));

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE restore_primal(SleqpRestoration* restoration,
                             const SleqpSparseVec* source,
                             SleqpSparseVec* target)
{
  int num_fixed_vars;
  int* fixed_var_indices;
  double* fixed_var_values;

  SLEQP_CALL(sleqp_preprocessing_state_fixed_variables(restoration->preprocessing_state,
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

static
SLEQP_RETCODE restore_converted_bound(SleqpConvertedBound* converted_bound,
                                      SLEQP_ACTIVE_STATE* var_state,
                                      SLEQP_ACTIVE_STATE* cons_state)
{
  const bool bound_flip = (converted_bound->factor < 0);

  switch (*var_state)
  {
  case SLEQP_INACTIVE:
    (*cons_state) = SLEQP_INACTIVE;
    break;
  case SLEQP_ACTIVE_BOTH:
    (*cons_state) = SLEQP_ACTIVE_BOTH;
    (*var_state) = SLEQP_ACTIVE_BOTH;
    break;
  case SLEQP_ACTIVE_UPPER:
    (*cons_state) = bound_flip ? SLEQP_ACTIVE_LOWER : SLEQP_ACTIVE_UPPER;
    (*var_state) = SLEQP_INACTIVE;
    break;
  case SLEQP_ACTIVE_LOWER:
    (*cons_state) = bound_flip ? SLEQP_ACTIVE_UPPER : SLEQP_ACTIVE_LOWER;
    (*var_state) = SLEQP_INACTIVE;
    break;
  }

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE store_working_set(const SleqpRestoration* restoration,
                                SleqpWorkingSet* working_set)
{
  SleqpProblem* problem = restoration->original_problem;

  const int num_variables = sleqp_problem_num_variables(problem);
  const int num_constraints = sleqp_problem_num_constraints(problem);

  for(int j = 0; j < num_variables; ++j)
  {
    SLEQP_ACTIVE_STATE state = restoration->working_var_states[j];

    if(state != SLEQP_INACTIVE)
    {
      SLEQP_CALL(sleqp_working_set_add_variable(working_set,
                                                j,
                                                state));
    }
  }

  for(int i = 0; i < num_constraints; ++i)
  {
    SLEQP_ACTIVE_STATE state = restoration->working_cons_states[i];

    if(state != SLEQP_INACTIVE)
    {
      SLEQP_CALL(sleqp_working_set_add_constraint(working_set,
                                                  i,
                                                  state));
    }
  }

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE reset_working_set_states(SleqpRestoration* restoration)
{
  SleqpProblem* problem = restoration->original_problem;

  const int num_variables = sleqp_problem_num_variables(problem);
  const int num_constraints = sleqp_problem_num_constraints(problem);

  for(int i = 0; i < num_constraints; ++i)
  {
    restoration->working_cons_states[i] = SLEQP_INACTIVE;
  }

  for(int j = 0; j < num_variables; ++j)
  {
    restoration->working_var_states[j] = SLEQP_INACTIVE;
  }

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE restore_working_set(SleqpRestoration* restoration,
                                  const SleqpWorkingSet* transformed,
                                  SleqpWorkingSet* original)
{
  SleqpPreprocessingState* preprocessing_state = restoration->preprocessing_state;
  SleqpProblem* problem = restoration->original_problem;

  const int num_variables = sleqp_problem_num_variables(problem);

  const int num_general = sleqp_problem_num_general_constraints(problem);
  const int num_linear = sleqp_problem_num_linear_constraints(problem);

  SleqpConvertedBound* converted_bounds;
  int num_converted_bounds;

  SLEQP_CALL(sleqp_preprocessing_state_converted_bounds(preprocessing_state,
                                                        &converted_bounds,
                                                        &num_converted_bounds));

  SLEQP_CALL(reset_working_set_states(restoration));

  SleqpVariableState* var_states = sleqp_preprocessing_state_variable_states(preprocessing_state);

  SleqpConstraintState* linear_cons_states = sleqp_preprocessing_state_linear_constraint_states(preprocessing_state);

  {

    int offset = 0;

    for(int j = 0; j < num_variables; ++j)
    {
      if(var_states[j].state == SLEQP_VAR_UNCHANGED)
      {
        restoration->working_var_states[j] = sleqp_working_set_get_variable_state(transformed,
                                                                                  j - offset);
      }
      else
      {
        assert(var_states[j].state == SLEQP_VAR_BOUNDFIXED);

        restoration->working_var_states[j] = SLEQP_ACTIVE_BOTH;

        ++offset;
      }
    }
  }

  {
    for(int i = 0; i < num_general; ++i)
    {
      restoration->working_cons_states[i] = sleqp_working_set_get_constraint_state(transformed,
                                                                                   i);
    }

    int offset = 0;

    for(int i = 0; i < num_linear; ++i)
    {
      const int i_general = i + num_general;

      if(linear_cons_states[i].state == SLEQP_CONS_UNCHANGED)
      {
        restoration->working_cons_states[i_general] = sleqp_working_set_get_constraint_state(transformed,
                                                                                             i_general - offset);
      }
      else if(linear_cons_states[i].state == SLEQP_CONS_REDUNDANT)
      {
        restoration->working_cons_states[i_general] = SLEQP_INACTIVE;
      }
      else
      {
        assert(linear_cons_states[i].state == SLEQP_CONS_BOUNDCONVERTED);

        SleqpConvertedBound* converted_bound = converted_bounds + linear_cons_states[i].bound;

        assert(converted_bound->constraint == i);

        SLEQP_CALL(restore_converted_bound(converted_bound,
                                           restoration->working_var_states + converted_bound->variable,
                                           restoration->working_cons_states + i_general));
      }
    }
  }

  store_working_set(restoration, original);

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE reset_duals(SleqpRestoration* restoration)
{
  SleqpProblem* problem = restoration->original_problem;

  const int num_variables = sleqp_problem_num_variables(problem);
  const int num_constraints = sleqp_problem_num_constraints(problem);

  for(int j = 0; j < num_variables; ++j)
  {
    restoration->var_dual[j] = SLEQP_INACTIVE;
  }

  for(int i = 0; i < num_constraints; ++i)
  {
    restoration->cons_dual[i] = SLEQP_INACTIVE;
  }

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE correct_dual_residuals(SleqpRestoration* restoration,
                                     SleqpIterate* original)
{
  SleqpPreprocessingState* preprocessing_state = restoration->preprocessing_state;

  SleqpProblem* problem = restoration->original_problem;

  const int num_variables = sleqp_problem_num_variables(problem);

  SleqpSparseVec* residuals = restoration->stationarity_residuals;

  const double zero_eps = sleqp_params_get(restoration->params,
                                           SLEQP_PARAM_ZERO_EPS);

  SLEQP_CALL(sleqp_iterate_stationarity_residuals(problem,
                                                  original,
                                                  restoration->cache,
                                                  residuals,
                                                  zero_eps));

  SleqpVariableState* var_states = sleqp_preprocessing_state_variable_states(preprocessing_state);

  for(int k = 0; k < residuals->nnz; ++k)
  {
    const int j = residuals->indices[k];
    const double v = residuals->data[k];

    if(var_states[j].state == SLEQP_VAR_BOUNDFIXED)
    {
      assert(restoration->var_dual[j] == 0.);
      restoration->var_dual[j] = -v;
    }
  }

  SLEQP_CALL(sleqp_sparse_vector_from_raw(sleqp_iterate_get_vars_dual(original),
                                          restoration->var_dual,
                                          num_variables,
                                          zero_eps));

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE restore_duals(SleqpRestoration* restoration,
                            const SleqpIterate* transformed,
                            SleqpIterate* original)
{
  SleqpPreprocessingState* preprocessing_state = restoration->preprocessing_state;

  SleqpProblem* problem = restoration->original_problem;

  const int num_variables = sleqp_problem_num_variables(problem);
  const int num_constraints = sleqp_problem_num_constraints(problem);
  const int num_general = sleqp_problem_num_general_constraints(problem);

  const double zero_eps = sleqp_params_get(restoration->params,
                                           SLEQP_PARAM_ZERO_EPS);

  SleqpConvertedBound* converted_bounds;
  int num_converted_bounds;

  SLEQP_CALL(sleqp_preprocessing_state_converted_bounds(preprocessing_state,
                                                        &converted_bounds,
                                                        &num_converted_bounds));

  SLEQP_CALL(reset_duals(restoration));

  {
    SleqpVariableState* var_states = sleqp_preprocessing_state_variable_states(preprocessing_state);

    SleqpSparseVec* var_dual = sleqp_iterate_get_vars_dual(transformed);

    int offset = 0;
    for(int k = 0; k < var_dual->nnz; ++k)
    {
      int j = var_dual->indices[k];
      double v = var_dual->data[k];

      if(var_states[j].state == SLEQP_VAR_UNCHANGED)
      {
        restoration->var_dual[j + offset] = v;
      }
      else
      {
        assert(var_states[j].state == SLEQP_VAR_BOUNDFIXED);
        ++offset;
      }
    }
  }

  {
    SleqpConstraintState* linear_cons_states = sleqp_preprocessing_state_linear_constraint_states(preprocessing_state);

    SleqpSparseVec* cons_dual = sleqp_iterate_get_cons_dual(transformed);

    int offset = 0;

    for(int k = 0; k < cons_dual->nnz; ++k)
    {
      int i = cons_dual->indices[k];
      double v = cons_dual->data[k];

      if(i < num_general)
      {
        restoration->var_dual[i] = v;
        continue;
      }

      const int i_general = i;
      i -= num_general;

      if(linear_cons_states[i].state == SLEQP_CONS_UNCHANGED)
      {
        restoration->cons_dual[i_general + offset] = v;
      }
      else
      {
        ++offset;
      }
    }

    for(int k = 0; k < num_converted_bounds; ++k)
    {
      SleqpConvertedBound* converted_bound = converted_bounds + linear_cons_states[k].bound;

      const int i_general = converted_bound->constraint + num_general;

      restoration->cons_dual[i_general] = restoration->var_dual[converted_bound->variable] / converted_bound->factor;
      restoration->var_dual[converted_bound->variable] = 0.;
    }
  }

  SLEQP_CALL(sleqp_sparse_vector_from_raw(sleqp_iterate_get_vars_dual(original),
                                          restoration->var_dual,
                                          num_variables,
                                          zero_eps));

  SLEQP_CALL(sleqp_sparse_vector_from_raw(sleqp_iterate_get_cons_dual(original),
                                          restoration->cons_dual,
                                          num_constraints,
                                          zero_eps));

  const int num_fixed_variables = sleqp_preprocessing_state_num_fixed_variables(preprocessing_state);

  if(num_fixed_variables > 0)
  {
    correct_dual_residuals(restoration, original);
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_restoration_restore_iterate(SleqpRestoration* restoration,
                                                const SleqpIterate* transformed,
                                                SleqpIterate* original)
{
  SleqpProblem* problem = restoration->original_problem;

  SLEQP_CALL(restore_primal(restoration,
                            sleqp_iterate_get_primal(transformed),
                            sleqp_iterate_get_primal(original)));

  SLEQP_CALL(sleqp_set_and_evaluate(problem,
                                    original,
                                    SLEQP_VALUE_REASON_NONE));

  SLEQP_CALL(restore_working_set(restoration,
                                 sleqp_iterate_get_working_set(transformed),
                                 sleqp_iterate_get_working_set(original)));

  SLEQP_CALL(restore_duals(restoration,
                           transformed,
                           original));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_restoration_capture(SleqpRestoration* restoration)
{
  ++restoration->refcount;

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE restoration_free(SleqpRestoration** star)
{
  SleqpRestoration* restoration = *star;

  SLEQP_CALL(sleqp_sparse_vector_free(&restoration->stationarity_residuals));

  sleqp_free(&restoration->cache);

  sleqp_free(&restoration->cons_dual);
  sleqp_free(&restoration->var_dual);

  sleqp_free(&restoration->working_cons_states);
  sleqp_free(&restoration->working_var_states);

  SLEQP_CALL(sleqp_params_release(&restoration->params));

  SLEQP_CALL(sleqp_problem_release(&restoration->original_problem));

  SLEQP_CALL(sleqp_preprocessing_state_release(&restoration->preprocessing_state));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_restoration_release(SleqpRestoration** star)
{
  SleqpRestoration* restoration = *star;

  if(!restoration)
  {
    return  SLEQP_OKAY;
  }

  if(--restoration->refcount == 0)
  {
    SLEQP_CALL(restoration_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
