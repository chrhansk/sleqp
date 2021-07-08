#include "preprocessor.h"

#include <stdlib.h>

#include "cmp.h"
#include "log.h"
#include "mem.h"
#include "util.h"

#include "preprocessor/fixed_var_func.h"
#include "preprocessor/preprocessing.h"

#include "preprocessor/preprocessing_state.h"
#include "preprocessor/transform.h"
#include "preprocessor/restore.h"

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
  SleqpRestoration* restoration;

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
      bound_state |= SLEQP_LOWER_BOUND;
      improved_lower = true;
    }
  }

  if(sleqp_is_finite(preprocessor->linear_ub[i]))
  {
    if(ub < preprocessor->var_ub[j])
    {
      bound_state |= SLEQP_UPPER_BOUND;
      improved_upper = true;
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
  const int num_linear = sleqp_problem_num_linear_constraints(problem);

  SleqpSparseMatrix* linear_coeffs = sleqp_problem_linear_coeffs(problem);

  assert(sleqp_sparse_matrix_get_num_rows(linear_coeffs) == num_linear);
  assert(sleqp_sparse_matrix_get_num_cols(linear_coeffs) == num_variables);

  double* linear_data = sleqp_sparse_matrix_get_data(linear_coeffs);
  int* linear_rows = sleqp_sparse_matrix_get_rows(linear_coeffs);
  int* linear_cols = sleqp_sparse_matrix_get_cols(linear_coeffs);

  const double inf = sleqp_infinity();

  for(int i = 0; i < num_linear; ++i)
  {
    preprocessor->linear_min[i] = 0.;
    preprocessor->linear_max[i] = 0.;
  }

  for(int col = 0; col < num_variables; ++col)
  {
    for(int k = linear_cols[col]; k < linear_cols[col + 1]; ++k)
    {
      const int row = linear_rows[k];
      const double value = linear_data[k];

      if(value == 0.)
      {
        continue;
      }

      double lb = preprocessor->var_lb[col];
      double ub = preprocessor->var_ub[col];

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
SLEQP_RETCODE add_forcing_constraint(SleqpPreprocessor* preprocessor,
                                     int i,
                                     SleqpBoundState bound_state)
{
  SLEQP_CALL(sleqp_preprocessing_state_add_forcing_constraint(preprocessor->preprocessing_state,
                                                              i,
                                                              bound_state,
                                                              preprocessor->var_lb,
                                                              preprocessor->var_ub));

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE check_for_constraint_infeasibility(SleqpPreprocessor* preprocessor)
{
  SleqpProblem* problem = preprocessor->original_problem;

  SleqpPreprocessingState* state = preprocessor->preprocessing_state;

  SleqpConstraintState* linear_cons_states = sleqp_preprocessing_state_linear_constraint_states(state);

  const int num_linear_constraints = sleqp_problem_num_linear_constraints(problem);

  for(int i = 0; i < num_linear_constraints; ++i)
  {
    if(linear_cons_states[i].state != SLEQP_CONS_UNCHANGED)
    {
      continue;
    }

    if(sleqp_is_finite(preprocessor->linear_lb[i]))
    {
      const double bound_slack = preprocessor->linear_max[i] - preprocessor->linear_lb[i];

      if(bound_slack < 0.)
      {
        preprocessor->infeasible = true;
      }
      else if(bound_slack == 0.)
      {
        add_forcing_constraint(preprocessor, i, SLEQP_LOWER_BOUND);
        // add forcing constraint
      }
    }

    if(sleqp_is_finite(preprocessor->linear_ub[i]))
    {
      const double bound_slack = preprocessor->linear_ub[i] - preprocessor->linear_min[i];

      if(bound_slack < 0.)
      {
        preprocessor->infeasible = true;
      }
      else if(bound_slack == 0.)
      {
        add_forcing_constraint(preprocessor, i, SLEQP_UPPER_BOUND);
        // add forcing constraint
      }
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
      SLEQP_CALL(sleqp_preprocessing_state_fix_variable_to_bounds(state,
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
      if(preprocessor->linear_lb[i] > 0. || preprocessor->linear_ub[i] < 0)
      {
        preprocessor->infeasible = true;
      }
      else
      {
        SLEQP_CALL(sleqp_preprocessing_state_remove_linear_constraint(preprocessor->preprocessing_state,
                                                                      i));
      }
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

  SLEQP_CALL(sleqp_preprocessing_state_flush(preprocessor->preprocessing_state));

  SLEQP_CALL(sleqp_transformation_create(&preprocessor->transformation,
                                         preprocessor->preprocessing_state,
                                         params));

  SLEQP_CALL(sleqp_transformation_create_transformed_problem(preprocessor->transformation,
                                                             &preprocessor->transformed_problem));

  SLEQP_CALL(sleqp_restoration_create(&preprocessor->restoration,
                                      preprocessor->preprocessing_state,
                                      params));

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

SLEQP_RETCODE sleqp_preprocessor_transform_primal(SleqpPreprocessor* preprocessor,
                                                  const SleqpSparseVec* source,
                                                  SleqpSparseVec* target)
{
  SLEQP_CALL(sleqp_transformation_convert_primal(preprocessor->transformation,
                                                 source,
                                                 target));
  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_preprocessor_restore_iterate(SleqpPreprocessor* preprocessor,
                                                 const SleqpIterate* transformed_iterate,
                                                 SleqpIterate* original_iterate)
{
  SLEQP_CALL(sleqp_restoration_restore_iterate(preprocessor->restoration,
                                               transformed_iterate,
                                               original_iterate));

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

  SLEQP_CALL(sleqp_restoration_release(&preprocessor->restoration));

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
