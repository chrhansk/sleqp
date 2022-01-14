#include "factorization_umfpack.h"

#include "fail.h"
#include <assert.h>

#include <umfpack.h>

#include "defs.h"
#include "error.h"
#include "log.h"
#include "mem.h"

typedef struct UmfpackData
{
  SleqpSparseMatrix* matrix;

  void* numeric_factorization;
  void* symbolic_factorization;

  double control[UMFPACK_CONTROL];
  double info[UMFPACK_INFO];

  double* solution;
  double* rhs;

  int current_size;

} UmfpackData;

static SLEQP_RETCODE
umfpack_get_error_string(int value, const char** message)
{
  switch (value)
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
    (*message) = "<Unknown>";
  }

  return SLEQP_OKAY;
}

#define UMFPACK_CALL(x)                                                        \
  do                                                                           \
  {                                                                            \
    int umfpack_status = (x);                                                  \
                                                                               \
    if (umfpack_status != UMFPACK_OK)                                          \
    {                                                                          \
      const char* umfpack_error_string;                                        \
      SLEQP_CALL(                                                              \
        umfpack_get_error_string(umfpack_status, &umfpack_error_string));      \
                                                                               \
      sleqp_raise(SLEQP_INTERNAL_ERROR,                                        \
                  "Caught Umfpack error <%d> (%s)",                            \
                  umfpack_status,                                              \
                  umfpack_error_string);                                       \
    }                                                                          \
  } while (0)

static SLEQP_RETCODE
umfpack_factorization_free_factorizations(UmfpackData* umfpack)
{
  if (umfpack->numeric_factorization)
  {
    umfpack_di_free_numeric(&umfpack->numeric_factorization);

    umfpack->numeric_factorization = NULL;
  }

  if (umfpack->symbolic_factorization)
  {
    umfpack_di_free_symbolic(&umfpack->symbolic_factorization);

    umfpack->symbolic_factorization = NULL;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
umfpack_factorization_set_matrix(void* factorization_data,
                                 SleqpSparseMatrix* matrix)
{
  UmfpackData* umfpack = (UmfpackData*)factorization_data;

  const int num_cols = sleqp_sparse_matrix_num_cols(matrix);
  const int num_rows = sleqp_sparse_matrix_num_rows(matrix);

  assert(num_cols == num_rows);

  if (umfpack->current_size < num_rows)
  {
    SLEQP_CALL(sleqp_realloc(&umfpack->rhs, num_rows));
    SLEQP_CALL(sleqp_realloc(&umfpack->solution, num_rows));

    umfpack->current_size = num_rows;
  }

  /*
  if(sleqp_log_level() >= SLEQP_LOG_DEBUG)
  {
    factorization->control[UMFPACK_PRL] = 2;
  }
  */

  SLEQP_CALL(umfpack_factorization_free_factorizations(umfpack));

  assert(sleqp_sparse_matrix_is_quadratic(matrix));

  umfpack->matrix = matrix;

  UMFPACK_CALL(umfpack_di_symbolic(sleqp_sparse_matrix_num_cols(matrix),
                                   sleqp_sparse_matrix_num_rows(matrix),
                                   sleqp_sparse_matrix_cols(matrix),
                                   sleqp_sparse_matrix_rows(matrix),
                                   sleqp_sparse_matrix_data(matrix),
                                   &umfpack->symbolic_factorization,
                                   umfpack->control,
                                   umfpack->info));

  UMFPACK_CALL(umfpack_di_numeric(sleqp_sparse_matrix_cols(matrix),
                                  sleqp_sparse_matrix_rows(matrix),
                                  sleqp_sparse_matrix_data(matrix),
                                  umfpack->symbolic_factorization,
                                  &umfpack->numeric_factorization,
                                  umfpack->control,
                                  umfpack->info));

  if (umfpack->symbolic_factorization)
  {
    umfpack_di_free_symbolic(&umfpack->symbolic_factorization);
  }

  umfpack->symbolic_factorization = NULL;

  for (int i = 0; i < num_rows; ++i)
  {
    umfpack->rhs[i] = 0.;
  }

  /*
  if(sleqp_log_level() >= SLEQP_LOG_DEBUG)
  {
    umfpack_di_report_info(umfpack->control,
                           umfpack->info);
  }
  */

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
set_cache(double* cache, const SleqpSparseVec* vec)
{
  for (int k = 0; k < vec->nnz; ++k)
  {
    cache[vec->indices[k]] = vec->data[k];
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
reset_cache(double* cache, const SleqpSparseVec* vec)
{
  for (int k = 0; k < vec->nnz; ++k)
  {
    cache[vec->indices[k]] = 0.;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
umfpack_factorization_solve(void* factorization_data, const SleqpSparseVec* rhs)
{
  UmfpackData* umfpack = (UmfpackData*)factorization_data;

  SleqpSparseMatrix* matrix = umfpack->matrix;

  assert(rhs->dim == sleqp_sparse_matrix_num_rows(matrix));

  SLEQP_CALL(set_cache(umfpack->rhs, rhs));

  UMFPACK_CALL(umfpack_di_solve(UMFPACK_A,
                                sleqp_sparse_matrix_cols(matrix),
                                sleqp_sparse_matrix_rows(matrix),
                                sleqp_sparse_matrix_data(matrix),
                                umfpack->solution,
                                umfpack->rhs,
                                umfpack->numeric_factorization,
                                umfpack->control,
                                umfpack->info));

  SLEQP_CALL(reset_cache(umfpack->rhs, rhs));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
umfpack_factorization_condition_estimate(void* factorization_data,
                                         double* condition_estimate)
{
  UmfpackData* umfpack = (UmfpackData*)factorization_data;

  *(condition_estimate) = 1. / (umfpack->info[UMFPACK_RCOND]);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
umfpack_factorization_solution(void* factorization_data,
                               SleqpSparseVec* sol,
                               int begin,
                               int end,
                               double zero_eps)
{
  UmfpackData* umfpack = (UmfpackData*)factorization_data;

  assert(begin <= end);

  SLEQP_CALL(sleqp_sparse_vector_from_raw(sol,
                                          umfpack->solution + begin,
                                          end - begin,
                                          zero_eps));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
umfpack_factorization_free(void** star)
{
  UmfpackData* umfpack = (UmfpackData*)(*star);

  if (!umfpack)
  {
    return SLEQP_OKAY;
  }

  sleqp_free(&umfpack->rhs);
  sleqp_free(&umfpack->solution);

  SLEQP_CALL(umfpack_factorization_free_factorizations(umfpack));

  sleqp_free(&umfpack);

  *star = NULL;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
umfpack_data_create(UmfpackData** star)
{
  SLEQP_CALL(sleqp_malloc(star));

  UmfpackData* umfpack_data = *star;

  *umfpack_data = (UmfpackData){0};

  umfpack_di_defaults(umfpack_data->control);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_factorization_umfpack_create(SleqpFactorization** star,
                                   SleqpParams* params)
{

  SleqpFactorizationCallbacks callbacks
    = {.set_matrix         = umfpack_factorization_set_matrix,
       .solve              = umfpack_factorization_solve,
       .solution           = umfpack_factorization_solution,
       .condition_estimate = umfpack_factorization_condition_estimate,
       .free               = umfpack_factorization_free};

  UmfpackData* umfpack_data;

  SLEQP_CALL(umfpack_data_create(&umfpack_data));

  SLEQP_CALL(sleqp_factorization_create(star,
                                        SLEQP_FACT_UMFPACK_NAME,
                                        SLEQP_FACT_UMFPACK_VERSION,
                                        params,
                                        &callbacks,
                                        (void*)umfpack_data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_factorization_create_default(SleqpFactorization** star,
                                   SleqpParams* params)
{
  SLEQP_CALL(sleqp_factorization_umfpack_create(star, params));

  return SLEQP_OKAY;
}
