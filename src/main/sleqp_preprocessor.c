#include "sleqp_preprocessor.h"

#include <stdlib.h>

#include "sleqp_cmp.h"
#include "sleqp_mem.h"

typedef struct Entry
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

  double* var_lb;
  double* var_ub;

  double* linear_lb;
  double* linear_ub;

  double* linear_min;
  double* linear_max;

  int num_redundant_linear_cons;
  int* redundant_linear_cons;

  SleqpSparseVec* transformed_var_lb;
  SleqpSparseVec* transformed_var_ub;

  SleqpSparseVec* transformed_linear_lb;
  SleqpSparseVec* transformed_linear_ub;

  SleqpSparseMatrix* transformed_linear_coeffs;

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
SLEQP_RETCODE
remove_redundant_linear_constraint(SleqpPreprocessor* preprocessor,
                                   int i)
{
  preprocessor->redundant_linear_cons[preprocessor->num_redundant_linear_cons++] = i;

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

  if(sleqp_is_finite(preprocessor->linear_ub[i]))
  {
    preprocessor->var_ub[j] = SLEQP_MIN(preprocessor->var_ub[j],
                                        ub);
  }

  if(sleqp_is_finite(preprocessor->linear_lb[i]))
  {
    preprocessor->var_lb[j] = SLEQP_MAX(preprocessor->var_lb[j],
                                        lb);
  }

  remove_redundant_linear_constraint(preprocessor, i);

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
SLEQP_RETCODE create_transformed_problem(SleqpPreprocessor* preprocessor)
{
  SleqpProblem* problem = preprocessor->original_problem;
  SleqpParams* params = preprocessor->params;

  const int num_variables = sleqp_problem_num_variables(problem);

  const double zero_eps = sleqp_params_get(params, SLEQP_PARAM_ZERO_EPS);

  sleqp_log_debug("Preprocessing removed %d linear constraints",
                  preprocessor->num_redundant_linear_cons);

  SLEQP_CALL(remove_vec_entries(sleqp_problem_linear_lb(problem),
                                preprocessor->transformed_linear_lb,
                                preprocessor->redundant_linear_cons,
                                preprocessor->num_redundant_linear_cons,
                                0));

  SLEQP_CALL(remove_vec_entries(sleqp_problem_linear_ub(problem),
                                preprocessor->transformed_linear_ub,
                                preprocessor->redundant_linear_cons,
                                preprocessor->num_redundant_linear_cons,
                                0));

  SLEQP_CALL(remove_matrix_rows(sleqp_problem_linear_coeffs(problem),
                                preprocessor->transformed_linear_coeffs,
                                preprocessor->redundant_linear_cons,
                                preprocessor->num_redundant_linear_cons,
                                0));

  SLEQP_CALL(sleqp_sparse_vector_from_raw(preprocessor->transformed_var_lb,
                                          preprocessor->var_lb,
                                          num_variables,
                                          zero_eps));

  SLEQP_CALL(sleqp_sparse_vector_from_raw(preprocessor->transformed_var_ub,
                                          preprocessor->var_ub,
                                          num_variables,
                                          zero_eps));

  SLEQP_CALL(sleqp_problem_create(&preprocessor->transformed_problem,
                                  sleqp_problem_func(problem),
                                  params,
                                  preprocessor->transformed_var_lb,
                                  preprocessor->transformed_var_ub,
                                  sleqp_problem_general_lb(problem),
                                  sleqp_problem_general_ub(problem),
                                  preprocessor->transformed_linear_coeffs,
                                  preprocessor->transformed_linear_lb,
                                  preprocessor->transformed_linear_ub));

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
SLEQP_RETCODE check_for_infeasibility(SleqpPreprocessor* preprocessor)
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

static int index_compare(const void* first,
                         const void* second)
{
  int i = *((int*) first);
  int j = *((int*) second);

  if(i < j)
  {
    return -1;
  }
  else if(i > j)
  {
    return 1;
  }

  return 0;
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

  qsort(preprocessor->redundant_linear_cons,
        preprocessor->num_redundant_linear_cons,
        sizeof(int),
        index_compare);

   SLEQP_CALL(check_for_infeasibility(preprocessor));

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

  SLEQP_CALL(sleqp_problem_capture(preprocessor->original_problem));

  SLEQP_CALL(sleqp_alloc_array(&preprocessor->linear_cons_counts, num_linear_constraints));
  SLEQP_CALL(sleqp_alloc_array(&preprocessor->linear_entries, num_linear_constraints));

  SLEQP_CALL(sleqp_alloc_array(&preprocessor->redundant_linear_cons, num_linear_constraints));

  SLEQP_CALL(sleqp_alloc_array(&preprocessor->var_lb, num_variables));
  SLEQP_CALL(sleqp_alloc_array(&preprocessor->var_ub, num_variables));

  SLEQP_CALL(sleqp_alloc_array(&preprocessor->linear_lb, num_linear_constraints));
  SLEQP_CALL(sleqp_alloc_array(&preprocessor->linear_ub, num_linear_constraints));

  SLEQP_CALL(sleqp_alloc_array(&preprocessor->linear_min, num_linear_constraints));
  SLEQP_CALL(sleqp_alloc_array(&preprocessor->linear_max, num_linear_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&preprocessor->transformed_var_lb,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&preprocessor->transformed_var_ub,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&preprocessor->transformed_linear_lb,
                                              num_linear_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&preprocessor->transformed_linear_ub,
                                              num_linear_constraints));

  SLEQP_CALL(sleqp_sparse_matrix_create(&preprocessor->transformed_linear_coeffs,
                                        num_linear_constraints,
                                        num_variables,
                                        0));

  SLEQP_CALL(transform_problem(preprocessor));

  return SLEQP_OKAY;
}

SLEQP_PREPROCESSING_RESULT sleqp_preprocessor_result(SleqpPreprocessor* preprocessor)
{
  if(preprocessor->infeasible)
  {
    return SLEQP_PREPROCESSING_RESULT_INFEASIBLE;
  }
  if(preprocessor->num_redundant_linear_cons > 0)
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
                                preprocessor->redundant_linear_cons,
                                preprocessor->num_redundant_linear_cons,
                                offset));

  SLEQP_CALL(remove_vec_entries(sleqp_iterate_get_cons_dual(iterate),
                                sleqp_iterate_get_cons_dual(transformed_iterate),
                                preprocessor->redundant_linear_cons,
                                preprocessor->num_redundant_linear_cons,
                                offset));

  SLEQP_CALL(remove_matrix_rows(sleqp_iterate_get_cons_jac(iterate),
                                sleqp_iterate_get_cons_jac(transformed_iterate),
                                preprocessor->redundant_linear_cons,
                                preprocessor->num_redundant_linear_cons,
                                offset));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_preprocessor_restore_iterate(SleqpPreprocessor* preprocessor,
                                                 const SleqpIterate* transformed_iterate,
                                                 SleqpIterate* original_iterate)
{
  return SLEQP_OKAY;
}

static SLEQP_RETCODE preprocessor_free(SleqpPreprocessor** star)
{
  SleqpPreprocessor* preprocessor = *star;

  SLEQP_CALL(sleqp_problem_release(&preprocessor->transformed_problem));

  SLEQP_CALL(sleqp_sparse_matrix_release(&preprocessor->transformed_linear_coeffs));

  SLEQP_CALL(sleqp_sparse_vector_free(&preprocessor->transformed_linear_ub));
  SLEQP_CALL(sleqp_sparse_vector_free(&preprocessor->transformed_linear_lb));

  SLEQP_CALL(sleqp_sparse_vector_free(&preprocessor->transformed_var_ub));
  SLEQP_CALL(sleqp_sparse_vector_free(&preprocessor->transformed_var_lb));

  sleqp_free(&preprocessor->linear_max);
  sleqp_free(&preprocessor->linear_min);

  sleqp_free(&preprocessor->linear_ub);
  sleqp_free(&preprocessor->linear_lb);

  sleqp_free(&preprocessor->var_ub);
  sleqp_free(&preprocessor->var_lb);

  sleqp_free(&preprocessor->redundant_linear_cons);

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
