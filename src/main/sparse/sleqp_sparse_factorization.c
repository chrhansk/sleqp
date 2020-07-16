#include "sleqp_sparse_factorization.h"

#include "assert.h"

#include <umfpack.h>

#include "sleqp_mem.h"

struct SleqpSparseFactorization
{
  SleqpSparseMatrix* matrix;

  void* numeric_factorization;
  void* symbolic_factorization;

  double control[UMFPACK_CONTROL];
  double info[UMFPACK_INFO];

  double* solution;
  double* rhs;
};

static SLEQP_RETCODE umfpack_get_error_string(int value,
                                              const char** message)
{
  switch(value)
  {
  case UMFPACK_ERROR_different_pattern:
    (*message) = "UMFPACK_ERROR_different_pattern";
    break;
  case UMFPACK_ERROR_invalid_Symbolic_object:
    (*message) = "UMFPACK_ERROR_invalid_Symbolic_object";
    break;
  case UMFPACK_WARNING_singular_matrix:
    (*message) = "UMFPACK_WARNING_singular_matrix";
    break;
  case UMFPACK_ERROR_n_nonpositive:
    (*message) = "UMFPACK_ERROR_n_nonpositive";
    break;
  case UMFPACK_ERROR_invalid_matrix:
    (*message) = "UMFPACK_ERROR_invalid_matrix";
    break;
  case UMFPACK_ERROR_argument_missing:
    (*message) = "MFPACK_ERROR_argument_missing";
    break;
  case UMFPACK_ERROR_out_of_memory:
    (*message) = "UMFPACK_ERROR_out_of_memory";
    break;
  case UMFPACK_ERROR_internal_error:
    (*message) = "UMFPACK_ERROR_internal_error";
    break;
  default:
    (*message) = "Unknown";
    return SLEQP_INTERNAL_ERROR;
  }

  return SLEQP_OKAY;
}


#define UMFPACK_CALL(x)                                                 \
  do                                                                    \
  {                                                                     \
    int umfpack_status = (x);                                           \
                                                                        \
    if(umfpack_status != UMFPACK_OK)                                    \
    {                                                                   \
      const char* umfpack_error_string;                                 \
      SLEQP_CALL(umfpack_get_error_string(umfpack_status,               \
                                          &umfpack_error_string));      \
                                                                        \
      sleqp_log_error("Caught Umfpack error <%d> (%s)",                 \
                      umfpack_status,                                   \
                      umfpack_error_string);                            \
                                                                        \
      switch(umfpack_status)                                            \
      {                                                                 \
      case UMFPACK_ERROR_invalid_matrix:                                \
      case UMFPACK_ERROR_argument_missing:                              \
        return SLEQP_ILLEGAL_ARGUMENT;                                  \
      case UMFPACK_ERROR_out_of_memory:                                 \
        return SLEQP_NOMEM;                                             \
      default:                                                          \
        return SLEQP_INTERNAL_ERROR;                                    \
      }                                                                 \
    }                                                                   \
  } while(0)

SLEQP_RETCODE sleqp_sparse_factorization_create(SleqpSparseFactorization** star,
                                                SleqpSparseMatrix* matrix)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpSparseFactorization* factorization = *star;

  *factorization = (SleqpSparseFactorization){0};

  umfpack_di_defaults(factorization->control);

  assert(sleqp_sparse_matrix_is_quadratic(matrix));

  factorization->matrix = matrix;

  UMFPACK_CALL(umfpack_di_symbolic(sleqp_sparse_matrix_get_num_cols(matrix),
                                   sleqp_sparse_matrix_get_num_rows(matrix),
                                   sleqp_sparse_matrix_get_cols(matrix),
                                   sleqp_sparse_matrix_get_rows(matrix),
                                   sleqp_sparse_matrix_get_data(matrix),
                                   &factorization->symbolic_factorization,
                                   factorization->control,
                                   factorization->info));

  UMFPACK_CALL(umfpack_di_numeric(sleqp_sparse_matrix_get_cols(matrix),
                                  sleqp_sparse_matrix_get_rows(matrix),
                                  sleqp_sparse_matrix_get_data(matrix),
                                  factorization->symbolic_factorization,
                                  &factorization->numeric_factorization,
                                  factorization->control,
                                  factorization->info));

  if(factorization->symbolic_factorization)
  {
    umfpack_di_free_symbolic(&factorization->symbolic_factorization);
  }

  factorization->symbolic_factorization = NULL;

  const int num_cols = sleqp_sparse_matrix_get_num_cols(matrix);
  const int num_rows = sleqp_sparse_matrix_get_num_rows(matrix);

  SLEQP_CALL(sleqp_calloc(&factorization->rhs, num_rows));

  for(int i = 0; i < num_rows; ++i)
  {
    factorization->rhs[i]  = 0.;
  }

  SLEQP_CALL(sleqp_calloc(&factorization->solution, num_cols));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE set_cache(double* cache,
                               SleqpSparseVec* vec)
{
  for(int k = 0; k < vec->nnz; ++k)
  {
    cache[vec->indices[k]] = vec->data[k];
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE reset_cache(double* cache,
                                 SleqpSparseVec* vec)
{
  for(int k = 0; k < vec->nnz; ++k)
  {
    cache[vec->indices[k]] = 0.;
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_factorization_solve(SleqpSparseFactorization* factorization,
                                               SleqpSparseVec* rhs)
{
  SleqpSparseMatrix* matrix = factorization->matrix;

  assert(rhs->dim == sleqp_sparse_matrix_get_num_rows(matrix));

  SLEQP_CALL(set_cache(factorization->rhs, rhs));

  UMFPACK_CALL(umfpack_di_solve(UMFPACK_A,
                                sleqp_sparse_matrix_get_cols(matrix),
                                sleqp_sparse_matrix_get_rows(matrix),
                                sleqp_sparse_matrix_get_data(matrix),
                                factorization->solution,
                                factorization->rhs,
                                factorization->numeric_factorization,
                                factorization->control,
                                factorization->info));

  SLEQP_CALL(reset_cache(factorization->rhs, rhs));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_factorization_get_condition_estimate(SleqpSparseFactorization* factorization,
                                                                double* condition_estimate)
{
  *(condition_estimate) = 1. / (factorization->info[UMFPACK_RCOND]);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_factorization_get_sol(SleqpSparseFactorization* factorization,
                                                 SleqpSparseVec* sol,
                                                 int begin,
                                                 int end,
                                                 double eps)
{
  assert(begin <= end);

  SLEQP_CALL(sleqp_sparse_vector_from_raw(sol,
                                          factorization->solution + begin,
                                          end - begin,
                                          eps));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_factorization_free(SleqpSparseFactorization** star)
{
  SleqpSparseFactorization* factorization = *star;

  if(!factorization)
  {
    return SLEQP_OKAY;
  }

  sleqp_free(&factorization->rhs);
  sleqp_free(&factorization->solution);

  if(factorization->numeric_factorization)
  {
    umfpack_di_free_numeric(&factorization->numeric_factorization);

    factorization->numeric_factorization = NULL;
  }

  if(factorization->symbolic_factorization)
  {
    umfpack_di_free_symbolic(&factorization->symbolic_factorization);

    factorization->symbolic_factorization = NULL;
  }

  sleqp_free(star);

  return SLEQP_OKAY;
}
