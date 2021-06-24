#include "sleqp_preprocessor.h"

#include <stdlib.h>

#include "sleqp_cmp.h"
#include "sleqp_mem.h"
#include "sleqp_util.h"

#include "preprocessor/sleqp_fixed_var_func.h"
#include "preprocessor/sleqp_preprocessing.h"

#include "preprocessor/sleqp_preprocessing_state.h"
#include "preprocessor/sleqp_transform.h"

typedef struct
{
  int row;
  int col;
  double value;
} Entry;


struct SleqpPreprocessor
{
  int refcount;

  SleqpParams* params;
  SleqpProblem* original_problem;

  int* linear_cons_counts;
  Entry* linear_entries;

  // dense version of variable bounds
  double* var_lb;
  double* var_ub;

  // dense version of linear bounds
  double* linear_lb;
  double* linear_ub;

  // actual bounds of linear constraints
  double* linear_min;
  double* linear_max;

  double* cons_dual_dense;
  SLEQP_ACTIVE_STATE* cons_state_dense;

  SleqpPreprocessingState* preprocessing_state;
  SleqpTransformation* transformation;

  int* removed_linear_cons;
  int num_removed_linear_cons;

  SleqpSparseVec* cache;

  bool infeasible;

  SleqpProblem* transformed_problem;
};

static
SLEQP_RETCODE compute_cons_counts(SleqpPreprocessor* preprocessor)
{
  SleqpProblem* problem = preprocessor->original_problem;

  const int num_variables = sleqp_problem_num_variables(problem);
  const int num_linear_constraints = sleqp_problem_num_linear_constraints(problem);

  SleqpSparseMatrix* linear_coeffs = sleqp_problem_linear_coeffs(problem);

  assert(sleqp_sparse_matrix_get_num_rows(linear_coeffs) == num_linear_constraints);
  assert(sleqp_sparse_matrix_get_num_cols(linear_coeffs) == num_variables);

  double* linear_data = sleqp_sparse_matrix_get_data(linear_coeffs);
  int* linear_rows = sleqp_sparse_matrix_get_rows(linear_coeffs);
  int* linear_cols = sleqp_sparse_matrix_get_cols(linear_coeffs);

  for(int i = 0; i < num_linear_constraints; ++i)
  {
    preprocessor->linear_cons_counts[i] = 0;
  }

  for(int col = 0; col < num_variables; ++col)
  {
    for(int k = linear_cols[col]; k < linear_cols[col + 1]; ++k)
    {
      const int row = linear_rows[k];
      const double value  = linear_data[k];

      if(linear_data[k] != 0.)
      {
        ++preprocessor->linear_cons_counts[linear_rows[k]];

        preprocessor->linear_entries[row] =
          (Entry) {.row = row, .col = col, .value = value};
      }
    }
  }

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE convert_linear_constraint_to_bound(SleqpPreprocessor* preprocessor,
                                                 int i)
{
  Entry* entry = preprocessor->linear_entries + i;

  double ub = preprocessor->linear_ub[i] / entry->value;
  double lb = preprocessor->linear_lb[i] / entry->value;

  assert(sleqp_is_finite(entry->value));

  if(entry->value < 0)
  {
    double t = ub;
    ub = lb;
    lb = t;
  }

  const int j = entry->col;

  bool improved_upper = false, improved_lower = false;

  SleqpBoundState bound_state = 0;

  if(sleqp_is_finite(preprocessor->linear_lb[i]))
  {
    if(lb > preprocessor->var_lb[j])
    {
      bound_state |= SleqpLowerBound;
      improved_lower = true;
      preprocessor->var_lb[j] = lb;
    }
  }

  if(sleqp_is_finite(preprocessor->linear_ub[i]))
  {
    if(ub < preprocessor->var_ub[j])
    {
      bound_state |= SleqpUpperBound;
      improved_upper = true;
      preprocessor->var_ub[j] = ub;
    }
  }

  const bool improved = (improved_lower || improved_upper);

  if(improved)
  {
    SLEQP_CALL(sleqp_preprocessing_state_convert_linear_constraint_to_bound(preprocessor->preprocessing_state,
                                                                            i,
                                                                            j,
                                                                            entry->value,
                                                                            lb,
                                                                            ub,
                                                                            bound_state));
  }
  else
  {
    SLEQP_CALL(sleqp_preprocessing_state_remove_linear_constraint(preprocessor->preprocessing_state,
                                                                  i));
  }

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE compute_linear_bounds(SleqpPreprocessor* preprocessor)
{
  SleqpProblem* problem = preprocessor->original_problem;

  const int num_variables = sleqp_problem_num_variables(problem);
  const int num_linear_constraints = sleqp_problem_num_linear_constraints(problem);

  SleqpSparseMatrix* linear_coeffs = sleqp_problem_linear_coeffs(problem);

  assert(sleqp_sparse_matrix_get_num_rows(linear_coeffs) == num_linear_constraints);
  assert(sleqp_sparse_matrix_get_num_cols(linear_coeffs) == num_variables);

  double* linear_data = sleqp_sparse_matrix_get_data(linear_coeffs);
  int* linear_rows = sleqp_sparse_matrix_get_rows(linear_coeffs);
  int* linear_cols = sleqp_sparse_matrix_get_cols(linear_coeffs);

  const double inf = sleqp_infinity();

  for(int i = 0; i < num_linear_constraints; ++i)
  {
    preprocessor->linear_min[i] = 0.;
    preprocessor->linear_max[i] = 0.;
  }

  for(int col = 0; col < num_variables; ++col)
  {
    for(int k = linear_cols[col]; k < linear_cols[col + 1]; ++k)
    {
      const int row = linear_rows[k];
      const double value  = linear_data[k];

      if(value == 0.)
      {
        continue;
      }

      double lb = preprocessor->var_lb[row];
      double ub = preprocessor->var_ub[row];

      if(value > 0.)
      {
        if(sleqp_is_finite(ub))
        {
          if(sleqp_is_finite(preprocessor->linear_max[row]))
          {
            preprocessor->linear_max[row] += value * ub;
          }
        }
        else
        {
          preprocessor->linear_max[row] = inf;
        }

        if(sleqp_is_finite(lb))
        {
          if(sleqp_is_finite(preprocessor->linear_min[row]))
          {
            preprocessor->linear_min[row] += value * lb;
          }
        }
        else
        {
          preprocessor->linear_min[row] = -inf;
        }
      }
      else if(value < 0.)
      {
        if(sleqp_is_finite(ub))
        {
          if(sleqp_is_finite(preprocessor->linear_min[row]))
          {
            preprocessor->linear_min[row] += value * ub;
          }
        }
        else
        {
          preprocessor->linear_min[row] = -inf;
        }

        if(sleqp_is_finite(ub))
        {
          if(sleqp_is_finite(preprocessor->linear_max[row]))
          {
            preprocessor->linear_max[row] += value*lb;
          }
        }
        else
        {
          preprocessor->linear_max[row] = inf;
        }
      }
    }
  }

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE check_for_constraint_infeasibility(SleqpPreprocessor* preprocessor)
{
  SleqpProblem* problem = preprocessor->original_problem;

  const int num_linear_constraints = sleqp_problem_num_linear_constraints(problem);

  for(int i = 0; i < num_linear_constraints; ++i)
  {
    if(!sleqp_is_infinite(-preprocessor->linear_lb[i]) &&
       preprocessor->linear_max[i] < preprocessor->linear_lb[i])
    {
      preprocessor->infeasible = true;
      break;
    }

    if(!sleqp_is_infinite(preprocessor->linear_ub[i]) &&
       preprocessor->linear_min[i] > preprocessor->linear_ub[i])
    {
      preprocessor->infeasible = true;
      break;
    }
  }

  return SLEQP_OKAY;
}


static
SLEQP_RETCODE fix_variables_by_bounds(SleqpPreprocessor* preprocessor)
{
  SleqpProblem* problem = preprocessor->original_problem;

  SleqpPreprocessingState* state = preprocessor->preprocessing_state;

  const int num_variables = sleqp_problem_num_variables(problem);

  for(int j = 0; j < num_variables; ++j)
  {
    if(preprocessor->var_lb[j] == preprocessor->var_ub[j])
    {
      SLEQP_CALL(sleqp_preprocessing_state_fix_variable(state,
                                                        j,
                                                        preprocessor->var_lb[j]));
    }
  }

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE remove_redundant_constraints(SleqpPreprocessor* preprocessor)
{
  SLEQP_CALL(compute_cons_counts(preprocessor));

  SleqpProblem* problem = preprocessor->original_problem;

  SLEQP_CALL(sleqp_sparse_vector_to_raw(sleqp_problem_linear_lb(problem),
                                        preprocessor->linear_lb));

  SLEQP_CALL(sleqp_sparse_vector_to_raw(sleqp_problem_linear_ub(problem),
                                        preprocessor->linear_ub));

  const int num_linear_constraints = sleqp_problem_num_linear_constraints(problem);

  for(int i = 0; i < num_linear_constraints; ++i)
  {
    const int count = preprocessor->linear_cons_counts[i];

    if(sleqp_is_infinite(-preprocessor->linear_lb[i]) &&
       sleqp_is_infinite(preprocessor->linear_ub[i]))
    {
      SLEQP_CALL(sleqp_preprocessing_state_remove_linear_constraint(preprocessor->preprocessing_state,
                                                                    i));
    }
    else if(count == 0)
    {
      SLEQP_CALL(sleqp_preprocessing_state_remove_linear_constraint(preprocessor->preprocessing_state,
                                                                    i));
    }
    else if(count == 1)
    {
      SLEQP_CALL(convert_linear_constraint_to_bound(preprocessor, i));
    }
  }

  SLEQP_CALL(compute_linear_bounds(preprocessor));

  for (int i = 0; i < num_linear_constraints; ++i)
  {
    const int count = preprocessor->linear_cons_counts[i];

    if(count <= 1)
    {
      continue;
    }

    if(!sleqp_is_infinite(-preprocessor->linear_min[i]) &&
       !sleqp_is_infinite(preprocessor->linear_max[i]))
    {
      if(preprocessor->linear_lb[i] <= preprocessor->linear_min[i] &&
         preprocessor->linear_ub[i] >= preprocessor->linear_max[i])
      {
        SLEQP_CALL(sleqp_preprocessing_state_remove_linear_constraint(preprocessor->preprocessing_state,
                                                                      i));
      }
    }
  }

  SLEQP_CALL(check_for_constraint_infeasibility(preprocessor));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_preprocessor_create(SleqpPreprocessor** star,
                                        SleqpProblem* problem,
                                        SleqpParams* params)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpPreprocessor* preprocessor = *star;

  *preprocessor = (SleqpPreprocessor) {0};

  preprocessor->refcount = 1;

  preprocessor->params = params;
  SLEQP_CALL(sleqp_params_capture(preprocessor->params));

  preprocessor->original_problem = problem;

  const int num_variables = sleqp_problem_num_variables(problem);
  const int num_linear_constraints = sleqp_problem_num_linear_constraints(problem);
  const int num_constraints = sleqp_problem_num_constraints(problem);

  SLEQP_CALL(sleqp_problem_capture(preprocessor->original_problem));

  SLEQP_CALL(sleqp_alloc_array(&preprocessor->linear_cons_counts, num_linear_constraints));
  SLEQP_CALL(sleqp_alloc_array(&preprocessor->linear_entries, num_linear_constraints));

  SLEQP_CALL(sleqp_preprocessing_state_create(&preprocessor->preprocessing_state,
                                              problem));

  SLEQP_CALL(sleqp_alloc_array(&preprocessor->var_lb, num_variables));
  SLEQP_CALL(sleqp_alloc_array(&preprocessor->var_ub, num_variables));

  SLEQP_CALL(sleqp_alloc_array(&preprocessor->linear_lb, num_linear_constraints));
  SLEQP_CALL(sleqp_alloc_array(&preprocessor->linear_ub, num_linear_constraints));

  SLEQP_CALL(sleqp_alloc_array(&preprocessor->linear_min, num_linear_constraints));
  SLEQP_CALL(sleqp_alloc_array(&preprocessor->linear_max, num_linear_constraints));

  SLEQP_CALL(sleqp_alloc_array(&preprocessor->cons_dual_dense, num_constraints));

  SLEQP_CALL(sleqp_alloc_array(&preprocessor->cons_state_dense, num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&preprocessor->cache,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_to_raw(sleqp_problem_var_lb(problem),
                                        preprocessor->var_lb));

  SLEQP_CALL(sleqp_sparse_vector_to_raw(sleqp_problem_var_ub(problem),
                                        preprocessor->var_ub));

  SLEQP_CALL(fix_variables_by_bounds(preprocessor));

  SLEQP_CALL(remove_redundant_constraints(preprocessor));

  SLEQP_CALL(sleqp_transformation_create(&preprocessor->transformation,
                                         preprocessor->preprocessing_state,
                                         params));

  SLEQP_CALL(sleqp_transformation_create_transformed_problem(preprocessor->transformation,
                                                             &preprocessor->transformed_problem));

  return SLEQP_OKAY;
}

SLEQP_PREPROCESSING_RESULT sleqp_preprocessor_result(SleqpPreprocessor* preprocessor)
{
  if(preprocessor->infeasible)
  {
    return SLEQP_PREPROCESSING_RESULT_INFEASIBLE;
  }

  SleqpPreprocessingState* state = preprocessor->preprocessing_state;

  const int num_fixed_vars = sleqp_preprocessing_state_num_fixed_variables(state);

  const int num_removed_cons = sleqp_preprocessing_state_num_removed_linear_constraints(state);

  if(num_fixed_vars > 0 || num_removed_cons > 0)
  {
    return SLEQP_PREPROCESSING_RESULT_SUCCESS;
  }

  return SLEQP_PREPROCESSING_RESULT_FAILURE;
}

SleqpProblem* sleqp_preprocessor_transformed_problem(SleqpPreprocessor* preprocessor)
{
  return preprocessor->transformed_problem;
}

SLEQP_RETCODE sleqp_preprocessor_transform_primal(const SleqpSparseVec* source,
                                                  SleqpSparseVec* target)
{
  return SLEQP_OKAY;
}
/*
static
SLEQP_RETCODE restore_duals(SleqpPreprocessor* preprocessor,
                            const SleqpIterate* transformed_iterate,
                            SleqpIterate* original_iterate)
{
  SleqpProblem* problem = preprocessor->original_problem;
  SleqpParams* params = preprocessor->params;

  const double zero_eps = sleqp_params_get(params, SLEQP_PARAM_ZERO_EPS);

  const int num_general = sleqp_problem_num_general_constraints(problem);

  const int num_constraints = sleqp_problem_num_constraints(problem);

  for(int i = 0; i < num_constraints; ++i)
  {
    preprocessor->cons_dual_dense[i] = 0.;
  }

  {
    const SleqpSparseVec* transformed_vars_dual = sleqp_iterate_get_vars_dual(transformed_iterate);
    SleqpSparseVec* original_vars_dual = sleqp_iterate_get_vars_dual(original_iterate);

    SLEQP_CALL(sleqp_sparse_vector_clear(original_vars_dual));

    SLEQP_CALL(sleqp_sparse_vector_reserve(original_vars_dual,
                                           transformed_vars_dual->nnz));

    for(int k = 0; k < transformed_vars_dual->nnz; ++k)
    {
      const int j = transformed_vars_dual->indices[k];
      const double v = transformed_vars_dual->data[k];

      if(preprocessor->bound_converted_vars[j] == SLEQP_NONE)
      {
        SLEQP_CALL(sleqp_sparse_vector_push(original_vars_dual,
                                            j,
                                            v));
      }
      else
      {
        const int bound_index = preprocessor->bound_converted_vars[j];

        ConvertedBound* bound = preprocessor->converted_bounds + bound_index;

        preprocessor->cons_dual_dense[num_general + bound->row] = v / bound->factor;
      }
    }
  }

  {
    const SleqpSparseVec* transformed_cons_dual = sleqp_iterate_get_cons_dual(transformed_iterate);
    SleqpSparseVec* original_cons_dual = sleqp_iterate_get_cons_dual(original_iterate);

    int offset = 0;

    for(int k = 0; k < transformed_cons_dual->nnz; ++k)
    {
      int i = transformed_cons_dual->indices[k];
      const double v = transformed_cons_dual->data[k];

      if(i < num_general)
      {
        preprocessor->cons_dual_dense[i] = v;
        continue;
      }

      int linear_index = i - num_general;

      while(preprocessor->linear_cons_states[linear_index + offset].state != ConsUnchanged)
      {
        ++offset;
      }

      preprocessor->cons_dual_dense[num_general + linear_index + offset] = v;
    }

    SLEQP_CALL(sleqp_sparse_vector_from_raw(original_cons_dual,
                                            preprocessor->cons_dual_dense,
                                            num_constraints,
                                            zero_eps));
  }

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE restore_working_set(SleqpPreprocessor* preprocessor,
                                  const SleqpWorkingSet* transformed_working_set,
                                  SleqpWorkingSet* original_working_set)
{
  SleqpProblem* problem = preprocessor->original_problem;
  SleqpProblem* transformed_problem = preprocessor->transformed_problem;

  const int num_variables = sleqp_problem_num_variables(problem);
  const int num_general = sleqp_problem_num_general_constraints(problem);

  const int num_constraints = sleqp_problem_num_constraints(problem);

  const int num_transformed_constraints = sleqp_problem_num_linear_constraints(transformed_problem);

  const int num_converted_bounds = preprocessor->stats.num_bound_converted_cons;

  SLEQP_CALL(sleqp_working_set_reset(original_working_set));

  // initially mark all converted constraints as inactive
  for(int k = 0; k < num_converted_bounds; ++k)
  {
    preprocessor->converted_bounds[k].var_state = SLEQP_INACTIVE;
  }

  for(int i = 0; i < num_constraints; ++i)
  {
    preprocessor->cons_state_dense[i] = SLEQP_INACTIVE;
  }

  // Set variables states minus the converted bounds
  for(int j = 0; j < num_variables; ++j)
  {
    SLEQP_ACTIVE_STATE var_state = sleqp_working_set_get_variable_state(transformed_working_set,
                                                                        j);

    if(var_state == SLEQP_INACTIVE)
    {
      continue;
    }

    const int bound_index = preprocessor->bound_converted_vars[j];

    if(bound_index != SLEQP_NONE)
    {
      // Variable bound determined by a previously converted linear constraint.
      // Store the state and continue without making the variable active

      ConvertedBound* bound = preprocessor->converted_bounds + bound_index;

      bound->var_state = var_state;

      SLEQP_ACTIVE_STATE cons_state = SLEQP_INACTIVE;

      BoundState bound_state = bound->bound_state;

      assert(bound->factor != 0.);

      const double pos = (bound->factor > 0.);

      if(bound_state & Lower)
      {
        cons_state |= pos ? SLEQP_ACTIVE_LOWER : SLEQP_ACTIVE_UPPER;
      }
      else if(bound_state & Upper)
      {
        cons_state |= pos ? SLEQP_ACTIVE_UPPER : SLEQP_ACTIVE_LOWER;
      }

      preprocessor->cons_state_dense[num_general + bound->row] = cons_state;
    }
    else
    {
      SLEQP_CALL(sleqp_working_set_add_variable(original_working_set,
                                                j,
                                                var_state));
    }
  }

  // Copy general constraint states
  for(int i = 0; i < num_general; ++i)
  {
    SLEQP_ACTIVE_STATE cons_state = sleqp_working_set_get_constraint_state(transformed_working_set,
                                                                           i);

    preprocessor->cons_state_dense[i] = cons_state;
  }

  int offset = 0;

  for(int i = num_general; i < num_transformed_constraints; ++i)
  {
    SLEQP_ACTIVE_STATE cons_state = sleqp_working_set_get_constraint_state(transformed_working_set,
                                                                           i);

    const int linear_index = i - num_general;

    while(preprocessor->linear_cons_states[linear_index + offset].state != ConsUnchanged)
    {
      ++offset;
    }

    preprocessor->cons_state_dense[num_general + linear_index + offset] = cons_state;
  }

  for(int i = 0; i < num_constraints; ++i)
  {
    SLEQP_ACTIVE_STATE cons_state = preprocessor->cons_state_dense[i];

    if(cons_state != SLEQP_INACTIVE)
    {
      SLEQP_CALL(sleqp_working_set_add_constraint(original_working_set,
                                                  i,
                                                  cons_state));
    }

  }

  return SLEQP_OKAY;
}
*/

SLEQP_RETCODE sleqp_preprocessor_restore_iterate(SleqpPreprocessor* preprocessor,
                                                 const SleqpIterate* transformed_iterate,
                                                 SleqpIterate* original_iterate)
{
  /*
  SLEQP_CALL(sleqp_sparse_vector_copy(sleqp_iterate_get_primal(transformed_iterate),
                                      sleqp_iterate_get_primal(original_iterate)));

  SLEQP_CALL(sleqp_set_and_evaluate(preprocessor->original_problem,
                                    original_iterate,
                                    SLEQP_VALUE_REASON_NONE));

  SLEQP_CALL(sleqp_sparse_vector_copy(sleqp_iterate_get_vars_dual(transformed_iterate),
                                      sleqp_iterate_get_vars_dual(original_iterate)));

  SLEQP_CALL(restore_working_set(preprocessor,
                                 sleqp_iterate_get_working_set(transformed_iterate),
                                 sleqp_iterate_get_working_set(original_iterate)));

  SLEQP_CALL(restore_duals(preprocessor,
                           transformed_iterate,
                           original_iterate));
  */


  return SLEQP_OKAY;
}

static SLEQP_RETCODE preprocessor_free(SleqpPreprocessor** star)
{
  SleqpPreprocessor* preprocessor = *star;

  SLEQP_CALL(sleqp_problem_release(&preprocessor->transformed_problem));

  SLEQP_CALL(sleqp_sparse_vector_free(&preprocessor->cache));

  sleqp_free(&preprocessor->removed_linear_cons);

  sleqp_free(&preprocessor->cons_state_dense);

  sleqp_free(&preprocessor->cons_dual_dense);

  SLEQP_CALL(sleqp_preprocessing_state_release(&preprocessor->preprocessing_state));

  SLEQP_CALL(sleqp_transformation_release(&preprocessor->transformation));

  sleqp_free(&preprocessor->linear_max);
  sleqp_free(&preprocessor->linear_min);

  sleqp_free(&preprocessor->linear_ub);
  sleqp_free(&preprocessor->linear_lb);

  sleqp_free(&preprocessor->var_ub);
  sleqp_free(&preprocessor->var_lb);

  sleqp_free(&preprocessor->removed_linear_cons);

  sleqp_free(&preprocessor->linear_entries);
  sleqp_free(&preprocessor->linear_cons_counts);

  SLEQP_CALL(sleqp_problem_release(&preprocessor->original_problem));

  SLEQP_CALL(sleqp_params_release(&preprocessor->params));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_preprocessor_capture(SleqpPreprocessor* preprocessor)
{
  ++preprocessor->refcount;
  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_preprocessor_release(SleqpPreprocessor** star)
{
  SleqpPreprocessor* preprocessor = *star;

  if(!preprocessor)
  {
    return SLEQP_OKAY;
  }

  if(--preprocessor->refcount == 0)
  {
    SLEQP_CALL(preprocessor_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
