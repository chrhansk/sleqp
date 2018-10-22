#include "sleqp_sparse_factorization.h"

#include "assert.h"

#include <umfpack.h>

#include "sleqp_mem.h"

struct SleqpSparseFactorization
{
  SleqpSparseMatrix* matrix;

  void* numeric_factorization;

  double* solution;
  double* rhs;
};

static SLEQP_RETCODE umfpack_get_error_string(int value,
                                              const char** message)
{
  switch(value)
  {
  case UMFPACK_ERROR_invalid_matrix:
    (*message) = "UMFPACK_ERROR_invalid_matrix";
    break;
  case UMFPACK_ERROR_argument_missing:
    (*message) = "UMFPACK_ERROR_invalid_matrix";
    break;
  case UMFPACK_ERROR_out_of_memory:
    (*message) = "UMFPACK_ERROR_out_of_memory";
    break;
  default:
    (*message) = "Unknown";
    return SLEQP_INTERNAL_ERROR;
  }

  return SLEQP_OKAY;
}


#define UMFPACK_CALL(x)                                            \
  do                                                               \
  {                                                                \
  int status = (x);                                                \
                                                                   \
  if(status != UMFPACK_OK)                                         \
  {                                                                \
    const char* umfpack_error_string;                              \
    SLEQP_CALL(umfpack_get_error_string(status,                    \
                                        &umfpack_error_string));   \
                                                                   \
    sleqp_log_error("Caught Umfpack error <%d> (%s)",              \
                    status,                                        \
                    umfpack_error_string);                         \
                                                                   \
    switch (status)                                                \
    {                                                              \
    case UMFPACK_ERROR_invalid_matrix:                             \
    case UMFPACK_ERROR_argument_missing:                           \
      return SLEQP_INVALID;                                        \
    case UMFPACK_ERROR_out_of_memory:                              \
      return SLEQP_NOMEM;                                          \
    default:                                                       \
      return SLEQP_INTERNAL_ERROR;                                 \
    }                                                              \
  }                                                                \
  } while(0)

SLEQP_RETCODE sleqp_sparse_factorization_create(SleqpSparseFactorization** star,
                                                SleqpSparseMatrix* matrix)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpSparseFactorization* factorization = *star;

  assert(matrix->num_cols == matrix->num_rows);

  factorization->matrix = matrix;

  void* symbolic_factorization;

  UMFPACK_CALL(umfpack_di_symbolic(matrix->num_cols,
                                   matrix->num_rows,
                                   matrix->cols,
                                   matrix->rows,
                                   matrix->data,
                                   &symbolic_factorization,
                                   NULL,
                                   NULL));

  UMFPACK_CALL(umfpack_di_numeric(matrix->cols,
                                  matrix->rows,
                                  matrix->data,
                                  symbolic_factorization,
                                  &factorization->numeric_factorization,
                                  NULL,
                                  NULL));

  umfpack_di_free_symbolic(&symbolic_factorization);

  SLEQP_CALL(sleqp_calloc(&factorization->rhs,
                          matrix->num_rows));

  for(int i = 0; i < matrix->num_rows; ++i)
  {
    factorization->rhs[i]  = 0.;
  }

  SLEQP_CALL(sleqp_calloc(&factorization->solution,
                          matrix->num_cols));

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

  assert(rhs->dim == matrix->num_rows);

  SLEQP_CALL(set_cache(factorization->rhs, rhs));

  UMFPACK_CALL(umfpack_di_solve(UMFPACK_A,
                                matrix->cols,
                                matrix->rows,
                                matrix->data,
                                factorization->solution,
                                factorization->rhs,
                                factorization->numeric_factorization,
                                NULL,
                                NULL));

  SLEQP_CALL(reset_cache(factorization->rhs, rhs));

  /*
  SLEQP_CALL(sleqp_sparse_vector_from_raw(sol,
                                          factorization->solution,
                                          matrix->num_cols));
  */


  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_factorization_get_sol(SleqpSparseFactorization* factorization,
                                                 SleqpSparseVec* sol,
                                                 int begin,
                                                 int end)
{
  assert(begin <= end);

  SLEQP_CALL(sleqp_sparse_vector_from_raw(sol,
                                          factorization->solution + begin,
                                          end - begin));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_factorization_free(SleqpSparseFactorization** star)
{
  SleqpSparseFactorization* factorization = *star;

  sleqp_free(&factorization->rhs);
  sleqp_free(&factorization->solution);

  umfpack_di_free_numeric(&factorization->numeric_factorization);

  sleqp_free(star);

  return SLEQP_OKAY;
}
