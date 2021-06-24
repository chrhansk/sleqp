#include "sleqp_preprocessor.h"

#include <stdlib.h>

#include "sleqp_cmp.h"
#include "sleqp_mem.h"
#include "sleqp_util.h"

#include "preprocessor/sleqp_fixed_var_func.h"
#include "preprocessor/sleqp_preprocessing.h"

typedef struct
{
  int row;
  int col;
  double value;
} Entry;

typedef enum
{
  Lower = (1 << 1),
  Upper = (1 << 2),
  Both = (Upper | Lower)
} BoundState;

typedef struct
{
  int row, col;
  double factor;
  BoundState bound_state;
  SLEQP_ACTIVE_STATE var_state;
} ConvertedBound;

typedef enum
{
  ConsUnchanged,
  ConsRedundant,
  ConsConvertedToBound
} _ConstraintState;

typedef struct
{
  _ConstraintState state;
  int bound;
} ConstraintState;

typedef enum
{
  VarUnchanged,
  VarFixedByBound
} _VariableState;

typedef struct
{
  _VariableState state;
  double value;
} VariableState;

typedef struct
{
  int num_redundant_cons;
  int num_bound_converted_cons;
  int num_fixed_variables;
} PreprocessingStats;


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

  ConstraintState* linear_cons_states;
  int* bound_converted_vars;

  VariableState* var_states;

  double* cons_dual_dense;
  SLEQP_ACTIVE_STATE* cons_state_dense;

  ConvertedBound* converted_bounds;

  int* removed_linear_cons;
  int num_removed_linear_cons;

  SleqpSparseVec* transformed_var_lb;
  SleqpSparseVec* transformed_var_ub;

  SleqpSparseVec* transformed_linear_lb;
  SleqpSparseVec* transformed_linear_ub;

  SleqpSparseVec* cache;

  SleqpSparseMatrix* transformed_linear_coeffs;

  bool infeasible;

  SleqpProblem* transformed_problem;

  SleqpFunc* transformed_func;

  int num_fixed_vars;
  int* fixed_var_indices;
  double* fixed_var_values;

  PreprocessingStats stats;
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
SLEQP_RETCODE
remove_redundant_linear_constraint(SleqpPreprocessor* preprocessor,
                                   int i)
{
  assert(preprocessor->linear_cons_states[i].state == ConsUnchanged);

  ++preprocessor->stats.num_redundant_cons;

  preprocessor->linear_cons_states[i].state = ConsRedundant;

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE
add_converted_bound(SleqpPreprocessor* preprocessor,
                    int row,
                    int col,
                    double factor,
                    bool improved_lower,
                    bool improved_upper)
{
  assert(preprocessor->linear_cons_states[row].state == ConsUnchanged);

  BoundState bound_state = 0;

  if(improved_lower)
  {
    bound_state |= Lower;
  }
  if(improved_upper)
  {
    bound_state |= Upper;
  }

  const int bound_index = preprocessor->stats.num_bound_converted_cons;

  ConvertedBound* converted_bound = preprocessor->converted_bounds + bound_index;

  *converted_bound = (ConvertedBound) { .row = row, .col = col, .factor = factor, .bound_state = bound_state};

  preprocessor->linear_cons_states[row] = (ConstraintState) {.state = ConsConvertedToBound, .bound = bound_index};

  preprocessor->bound_converted_vars[col] = bound_index;

  ++preprocessor->stats.num_bound_converted_cons;

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE
convert_linear_constraint_to_bound(SleqpPreprocessor* preprocessor,
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

  if(sleqp_is_finite(preprocessor->linear_ub[i]))
  {
    if(ub < preprocessor->var_ub[j])
    {
      improved_upper = true;
      preprocessor->var_ub[j] = ub;
    }
  }

  if(sleqp_is_finite(preprocessor->linear_lb[i]))
  {
    if(lb > preprocessor->var_lb[j])
    {
      improved_lower = true;
      preprocessor->var_lb[j] = lb;
    }
  }

  const bool improved = (improved_lower || improved_upper);

  if(improved)
  {
    SLEQP_CALL(add_converted_bound(preprocessor,
                                   i,
                                   j,
                                   entry->value,
                                   improved_lower,
                                   improved_upper));
  }
  else
  {
    remove_redundant_linear_constraint(preprocessor, i);
  }

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE remove_vec_entries(const SleqpSparseVec* original,
                                 SleqpSparseVec* transformed,
                                 const int *entries,
                                 int num_entries,
                                 int offset)
{
  assert(num_entries <= original->dim);

  SLEQP_CALL(sleqp_sparse_vector_clear(transformed));

  SLEQP_CALL(sleqp_sparse_vector_resize(transformed,
                                        original->dim - num_entries));

  SLEQP_CALL(sleqp_sparse_vector_reserve(transformed,
                                         original->nnz));

  int i = offset;

  for(int k = 0; k < original->nnz; ++k)
  {
    if(i < num_entries && original->indices[k] == entries[i])
    {
      ++i;
    }
    else
    {
      SLEQP_CALL(sleqp_sparse_vector_push(transformed,
                                          original->indices[k] - i,
                                          original->data[k]));
    }
  }

  assert(i <= num_entries);

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE remove_matrix_rows(SleqpSparseMatrix* original,
                                 SleqpSparseMatrix* transformed,
                                 const int* entries,
                                 int num_entries,
                                 int offset)
{
  const int num_cols = sleqp_sparse_matrix_get_num_cols(original);
  const int num_rows = sleqp_sparse_matrix_get_num_rows(original);

  assert(num_entries <= num_rows);

  SLEQP_CALL(sleqp_sparse_matrix_clear(transformed));

  SLEQP_CALL(sleqp_sparse_matrix_resize(transformed,
                                        num_rows - num_entries,
                                        num_cols));

  SLEQP_CALL(sleqp_sparse_matrix_reserve(transformed,
                                         sleqp_sparse_matrix_get_nnz(original)));

  double* original_data = sleqp_sparse_matrix_get_data(original);
  int* original_rows = sleqp_sparse_matrix_get_rows(original);
  int* original_cols = sleqp_sparse_matrix_get_cols(original);

  for(int col = 0; col < num_cols; ++col)
  {
    int i = offset;

    SLEQP_CALL(sleqp_sparse_matrix_push_column(transformed,
                                               col));

    for(int k = original_cols[col]; k < original_cols[col + 1]; ++k)
    {
      const int row = original_rows[k];

      if(i < num_entries && row == entries[i])
      {
        ++i;
      }
      else
      {
        const double value = original_data[k];

        SLEQP_CALL(sleqp_sparse_matrix_push(transformed,
                                            row - i,
                                            col,
                                            value));
      }
    }

    assert(i <= num_entries);
  }

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE create_fixed_vars(SleqpPreprocessor* preprocessor)
{
  SleqpProblem* problem = preprocessor->original_problem;

  const int num_variables = sleqp_problem_num_variables(problem);

  const int num_fixed_vars = preprocessor->stats.num_fixed_variables;

  sleqp_free(&preprocessor->fixed_var_indices);
  sleqp_free(&preprocessor->fixed_var_values);

  SLEQP_CALL(sleqp_alloc_array(&preprocessor->fixed_var_indices, num_fixed_vars));
  SLEQP_CALL(sleqp_alloc_array(&preprocessor->fixed_var_values, num_fixed_vars));

  int i = 0;

  for(int j = 0; j < num_variables; ++j)
  {
    if(preprocessor->var_states[j].state != VarUnchanged)
    {
      assert(preprocessor->var_states[j].state == VarFixedByBound);

      preprocessor->fixed_var_indices[i] = j;
      preprocessor->fixed_var_values[i] = j;
      ++i;
    }
  }

  assert(i == num_fixed_vars);

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE create_transformed_func(SleqpPreprocessor* preprocessor,
                                      SleqpFunc** star)
{
  SleqpProblem* problem = preprocessor->original_problem;

  const int num_fixed_vars = preprocessor->stats.num_fixed_variables;

  SLEQP_CALL(sleqp_fixed_var_func_create(star,
                                         sleqp_problem_func(problem),
                                         num_fixed_vars,
                                         preprocessor->fixed_var_indices,
                                         preprocessor->fixed_var_values));

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE transform_var_bounds(SleqpPreprocessor* preprocessor)
{
  SleqpProblem* problem = preprocessor->original_problem;

  const int num_fixed_vars = preprocessor->stats.num_fixed_variables;

  const int num_variables = sleqp_problem_num_variables(problem);

  const bool has_fixed_vars = num_fixed_vars > 0;

  const double zero_eps = sleqp_params_get(preprocessor->params,
                                           SLEQP_PARAM_ZERO_EPS);

  if(has_fixed_vars)
  {
    SLEQP_CALL(sleqp_sparse_vector_resize(preprocessor->transformed_var_lb,
                                          num_variables - num_fixed_vars));

    SLEQP_CALL(sleqp_sparse_vector_from_raw(preprocessor->cache,
                                            preprocessor->var_lb,
                                            num_variables,
                                            zero_eps));

    SLEQP_CALL(sleqp_preprocessing_remove_fixed_entries(preprocessor->cache,
                                                        preprocessor->transformed_var_lb,
                                                        num_fixed_vars,
                                                        preprocessor->fixed_var_indices));

    SLEQP_CALL(sleqp_sparse_vector_resize(preprocessor->transformed_var_ub,
                                          num_variables - num_fixed_vars));

    SLEQP_CALL(sleqp_sparse_vector_from_raw(preprocessor->cache,
                                            preprocessor->var_ub,
                                            num_variables,
                                            zero_eps));

    SLEQP_CALL(sleqp_preprocessing_remove_fixed_entries(preprocessor->cache,
                                                        preprocessor->transformed_var_ub,
                                                        num_fixed_vars,
                                                        preprocessor->fixed_var_indices));

  }
  else
  {
    SLEQP_CALL(sleqp_sparse_vector_from_raw(preprocessor->transformed_var_lb,
                                            preprocessor->var_lb,
                                            num_variables,
                                            zero_eps));

    SLEQP_CALL(sleqp_sparse_vector_from_raw(preprocessor->transformed_var_ub,
                                            preprocessor->var_ub,
                                            num_variables,
                                            zero_eps));
  }


  return SLEQP_OKAY;
}

static
SLEQP_RETCODE create_transformed_problem(SleqpPreprocessor* preprocessor)
{
  SleqpProblem* problem = preprocessor->original_problem;
  SleqpParams* params = preprocessor->params;

  const bool has_fixed_vars = preprocessor->stats.num_fixed_variables > 0;
  SleqpFunc* transformed_func = NULL;

  SLEQP_CALL(create_fixed_vars(preprocessor));

  SLEQP_CALL(remove_vec_entries(sleqp_problem_linear_lb(problem),
                                preprocessor->transformed_linear_lb,
                                preprocessor->removed_linear_cons,
                                preprocessor->num_removed_linear_cons,
                                0));

  SLEQP_CALL(remove_vec_entries(sleqp_problem_linear_ub(problem),
                                preprocessor->transformed_linear_ub,
                                preprocessor->removed_linear_cons,
                                preprocessor->num_removed_linear_cons,
                                0));

  SLEQP_CALL(remove_matrix_rows(sleqp_problem_linear_coeffs(problem),
                                preprocessor->transformed_linear_coeffs,
                                preprocessor->removed_linear_cons,
                                preprocessor->num_removed_linear_cons,
                                0));

  SLEQP_CALL(transform_var_bounds(preprocessor));

  if(has_fixed_vars)
  {
    SLEQP_CALL(create_transformed_func(preprocessor, &transformed_func));
  }
  else
  {
    transformed_func = sleqp_problem_func(problem);
  }


  SLEQP_CALL(sleqp_problem_create(&preprocessor->transformed_problem,
                                  transformed_func,
                                  params,
                                  preprocessor->transformed_var_lb,
                                  preprocessor->transformed_var_ub,
                                  sleqp_problem_general_lb(problem),
                                  sleqp_problem_general_ub(problem),
                                  preprocessor->transformed_linear_coeffs,
                                  preprocessor->transformed_linear_lb,
                                  preprocessor->transformed_linear_ub));

  if(has_fixed_vars)
  {
    SLEQP_CALL(sleqp_func_release(&transformed_func));
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

/*
static
SLEQP_RETCODE fix_variables_by_bounds(SleqpPreprocessor* preprocessor)
{
  SleqpProblem* problem = preprocessor->original_problem;

  const int num_variables = sleqp_problem_num_variables(problem);

  for(int j = 0; j < num_variables; ++j)
  {
    if(preprocessor->var_lb[j] == preprocessor->var_ub[j])
    {
      preprocessor->var_states[j] = (VariableState) {.state = VarFixedByBound, .value = preprocessor->var_lb[j]};
      ++preprocessor->stats.num_fixed_variables;
    }
  }

  return SLEQP_OKAY;
}
*/

static
SLEQP_RETCODE collect_removed_linear_cons(SleqpPreprocessor* preprocessor)
{
  SleqpProblem* problem = preprocessor->original_problem;

  const int num_linear_constraints = sleqp_problem_num_linear_constraints(problem);

  const int num_removed_cons = preprocessor->stats.num_bound_converted_cons +
    preprocessor->stats.num_redundant_cons;

  preprocessor->num_removed_linear_cons = 0;

  SLEQP_CALL(sleqp_realloc(&preprocessor->removed_linear_cons,
                           num_removed_cons));

  for(int i = 0; i < num_linear_constraints; ++i)
  {
    if(preprocessor->linear_cons_states[i].state != ConsUnchanged)
    {
      preprocessor->removed_linear_cons[preprocessor->num_removed_linear_cons++] = i;
    }
  }

  assert(num_removed_cons == preprocessor->num_removed_linear_cons);

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE transform_problem(SleqpPreprocessor* preprocessor)
{
  SLEQP_CALL(compute_cons_counts(preprocessor));

  SleqpProblem* problem = preprocessor->original_problem;

  SLEQP_CALL(sleqp_sparse_vector_to_raw(sleqp_problem_var_lb(problem),
                                        preprocessor->var_lb));

  SLEQP_CALL(sleqp_sparse_vector_to_raw(sleqp_problem_var_ub(problem),
                                        preprocessor->var_ub));

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
      SLEQP_CALL(remove_redundant_linear_constraint(preprocessor, i));
    }
    else if(count == 0)
    {
      SLEQP_CALL(remove_redundant_linear_constraint(preprocessor, i));
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
        SLEQP_CALL(remove_redundant_linear_constraint(preprocessor, i));
      }
    }
  }

  SLEQP_CALL(check_for_constraint_infeasibility(preprocessor));

  SLEQP_CALL(collect_removed_linear_cons(preprocessor));

  // SLEQP_CALL(fix_variables_by_bounds(preprocessor));

  SLEQP_CALL(create_transformed_problem(preprocessor));

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

  SLEQP_CALL(sleqp_alloc_array(&preprocessor->linear_cons_states, num_linear_constraints));

  for(int i = 0; i < num_linear_constraints; ++i)
  {
    preprocessor->linear_cons_states[i] = (ConstraintState) {.state = ConsUnchanged, .bound = -1 };
  }

  SLEQP_CALL(sleqp_alloc_array(&preprocessor->converted_bounds,
                               num_linear_constraints));

  SLEQP_CALL(sleqp_alloc_array(&preprocessor->var_states, num_variables));

  for(int j = 0; j < num_variables; ++j)
  {
    preprocessor->var_states[j] = (VariableState) {.state = VarUnchanged};
  }

  SLEQP_CALL(sleqp_alloc_array(&preprocessor->var_lb, num_variables));
  SLEQP_CALL(sleqp_alloc_array(&preprocessor->var_ub, num_variables));

  SLEQP_CALL(sleqp_alloc_array(&preprocessor->linear_lb, num_linear_constraints));
  SLEQP_CALL(sleqp_alloc_array(&preprocessor->linear_ub, num_linear_constraints));

  SLEQP_CALL(sleqp_alloc_array(&preprocessor->linear_min, num_linear_constraints));
  SLEQP_CALL(sleqp_alloc_array(&preprocessor->linear_max, num_linear_constraints));

  SLEQP_CALL(sleqp_alloc_array(&preprocessor->bound_converted_vars, num_variables));

  for(int j = 0; j < num_variables; ++j)
  {
    preprocessor->bound_converted_vars[j] = SLEQP_NONE;
  }

  SLEQP_CALL(sleqp_alloc_array(&preprocessor->cons_dual_dense, num_constraints));

  SLEQP_CALL(sleqp_alloc_array(&preprocessor->cons_state_dense, num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&preprocessor->transformed_var_lb,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&preprocessor->transformed_var_ub,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&preprocessor->transformed_linear_lb,
                                              num_linear_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&preprocessor->transformed_linear_ub,
                                              num_linear_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&preprocessor->cache,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_matrix_create(&preprocessor->transformed_linear_coeffs,
                                        num_linear_constraints,
                                        num_variables,
                                        0));

  preprocessor->stats = (PreprocessingStats) {0};

  SLEQP_CALL(transform_problem(preprocessor));

  return SLEQP_OKAY;
}

SLEQP_PREPROCESSING_RESULT sleqp_preprocessor_result(SleqpPreprocessor* preprocessor)
{
  if(preprocessor->infeasible)
  {
    return SLEQP_PREPROCESSING_RESULT_INFEASIBLE;
  }

  const int num_removed_cons = preprocessor->stats.num_bound_converted_cons +
    preprocessor->stats.num_redundant_cons;

  if(num_removed_cons > 0)
  {
    return SLEQP_PREPROCESSING_RESULT_SUCCESS;
  }

  return SLEQP_PREPROCESSING_RESULT_FAILURE;
}

SleqpProblem* sleqp_preprocessor_transformed_problem(SleqpPreprocessor* preprocessor)
{
  return preprocessor->transformed_problem;
}

SLEQP_RETCODE sleqp_preprocessor_transform_iterate(SleqpPreprocessor* preprocessor,
                                                   const SleqpIterate* iterate,
                                                   SleqpIterate* transformed_iterate)
{
  SleqpProblem* problem = preprocessor->original_problem;

  const int offset = sleqp_problem_num_general_constraints(problem);

  SLEQP_CALL(sleqp_sparse_vector_copy(sleqp_iterate_get_primal(iterate),
                                      sleqp_iterate_get_primal(transformed_iterate)));

  SLEQP_CALL(remove_vec_entries(sleqp_iterate_get_cons_val(iterate),
                                sleqp_iterate_get_cons_val(transformed_iterate),
                                preprocessor->removed_linear_cons,
                                preprocessor->num_removed_linear_cons,
                                offset));

  SLEQP_CALL(remove_vec_entries(sleqp_iterate_get_cons_dual(iterate),
                                sleqp_iterate_get_cons_dual(transformed_iterate),
                                preprocessor->removed_linear_cons,
                                preprocessor->num_removed_linear_cons,
                                offset));

  SLEQP_CALL(remove_matrix_rows(sleqp_iterate_get_cons_jac(iterate),
                                sleqp_iterate_get_cons_jac(transformed_iterate),
                                preprocessor->removed_linear_cons,
                                preprocessor->num_removed_linear_cons,
                                offset));

  return SLEQP_OKAY;
}

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

SLEQP_RETCODE sleqp_preprocessor_restore_iterate(SleqpPreprocessor* preprocessor,
                                                 const SleqpIterate* transformed_iterate,
                                                 SleqpIterate* original_iterate)
{
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


  return SLEQP_OKAY;
}

static SLEQP_RETCODE preprocessor_free(SleqpPreprocessor** star)
{
  SleqpPreprocessor* preprocessor = *star;

  sleqp_free(&preprocessor->fixed_var_indices);
  sleqp_free(&preprocessor->fixed_var_values);

  SLEQP_CALL(sleqp_problem_release(&preprocessor->transformed_problem));

  SLEQP_CALL(sleqp_sparse_matrix_release(&preprocessor->transformed_linear_coeffs));

  SLEQP_CALL(sleqp_sparse_vector_free(&preprocessor->cache));

  SLEQP_CALL(sleqp_sparse_vector_free(&preprocessor->transformed_linear_ub));
  SLEQP_CALL(sleqp_sparse_vector_free(&preprocessor->transformed_linear_lb));

  SLEQP_CALL(sleqp_sparse_vector_free(&preprocessor->transformed_var_ub));
  SLEQP_CALL(sleqp_sparse_vector_free(&preprocessor->transformed_var_lb));

  sleqp_free(&preprocessor->removed_linear_cons);

  sleqp_free(&preprocessor->converted_bounds);

  sleqp_free(&preprocessor->cons_state_dense);

  sleqp_free(&preprocessor->cons_dual_dense);

  sleqp_free(&preprocessor->bound_converted_vars);

  sleqp_free(&preprocessor->linear_max);
  sleqp_free(&preprocessor->linear_min);

  sleqp_free(&preprocessor->linear_ub);
  sleqp_free(&preprocessor->linear_lb);

  sleqp_free(&preprocessor->var_ub);
  sleqp_free(&preprocessor->var_lb);

  sleqp_free(&preprocessor->removed_linear_cons);

  sleqp_free(&preprocessor->linear_cons_states);

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
