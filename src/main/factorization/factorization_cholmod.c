#include "factorization_cholmod.h"

#include "factorization/factorization.h"
#include "fail.h"
#include <assert.h>

#include <cholmod.h>
#include <cholmod_core.h>

#include "defs.h"
#include "error.h"
#include "log.h"
#include "mem.h"
#include "pub_log.h"
#include "pub_types.h"
#include "sparse/pub_sparse_matrix.h"

#define CHOLMOD_FLAGS (SLEQP_FACTORIZATION_PSD | SLEQP_FACTORIZATION_LOWER)

typedef struct SPQRData
{
  SleqpSparseMatrix* matrix;

  cholmod_common common;
  cholmod_sparse* sparse;

  cholmod_factor* factor;

  cholmod_dense* sol;
  cholmod_dense* rhs;

  int num_rows;
  int num_cols;

} CHOLMODData;

static void
report_error(int status, const char* file, int line, const char* message)
{
  if (status < 0)
  {
    sleqp_log_error("CHOLMOD: '%s'", message);
  }
  else
  {
    sleqp_log_warn("CHOLMOD: '%s'", message);
  }
}

static SLEQP_RETCODE
cholmod_error_string(int status, const char** message)
{
  switch (status)
  {
  case CHOLMOD_NOT_INSTALLED:
    (*message) = "method not installed";
  case CHOLMOD_OUT_OF_MEMORY:
    (*message) = "out of memory";
  case CHOLMOD_TOO_LARGE:
    (*message) = "integer overflow occured";
  case CHOLMOD_INVALID:
    (*message) = "invalid input";
  case CHOLMOD_GPU_PROBLEM:
    (*message) = "GPU fatal error";
  default:
    (*message) = "unknown error";
  }

  return SLEQP_OKAY;
}

#define CHOLMOD_ERROR_CHECK(common)                                            \
  do                                                                           \
  {                                                                            \
    if ((common)->status < 0)                                                  \
    {                                                                          \
      const char* message;                                                     \
      SLEQP_CALL(cholmod_error_string((common)->status, &message));            \
      sleqp_raise(SLEQP_INTERNAL_ERROR,                                        \
                  "Encountered error in SPQR: %s",                             \
                  message);                                                    \
    }                                                                          \
  } while (false)

static SLEQP_RETCODE
update_shape(CHOLMODData* cholmod_data, int num_rows, int num_cols, int nnz_max)
{
  cholmod_common* common = &(cholmod_data->common);

  if (num_rows != cholmod_data->num_rows)
  {
    cholmod_l_free_dense(&cholmod_data->rhs, common);

    cholmod_data->rhs
      = cholmod_l_allocate_dense(num_rows, 1, num_rows, CHOLMOD_REAL, common);

    CHOLMOD_ERROR_CHECK(common);

    for (int i = 0; i < num_rows; ++i)
    {
      ((double*)cholmod_data->rhs->x)[i] = 0.;
    }
  }

  if ((!cholmod_data->sparse) || (num_cols != cholmod_data->num_cols)
      || (num_rows != cholmod_data->num_rows))
  {
    cholmod_l_free_sparse(&cholmod_data->sparse, common);

    cholmod_data->sparse
      = cholmod_l_allocate_sparse(num_rows,
                                  num_cols,
                                  nnz_max,
                                  true, // sorted
                                  true, // packed
                                  -1,   // only lower part is stored
                                  CHOLMOD_REAL,
                                  common);

    CHOLMOD_ERROR_CHECK(common);
  }
  else
  {
    if (nnz_max > cholmod_data->sparse->nzmax)
    {
      cholmod_l_reallocate_sparse(nnz_max, cholmod_data->sparse, common);

      CHOLMOD_ERROR_CHECK(common);
    }
  }

  cholmod_data->num_cols = num_cols;
  cholmod_data->num_rows = num_rows;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
cholmod_factorization_set_matrix(void* factorization_data,
                                 SleqpSparseMatrix* matrix)
{
  CHOLMODData* cholmod_data = (CHOLMODData*)factorization_data;

  cholmod_common* common = &(cholmod_data->common);

  const int num_rows = sleqp_sparse_matrix_num_rows(matrix);
  const int num_cols = sleqp_sparse_matrix_num_cols(matrix);
  const int nnz      = sleqp_sparse_matrix_nnz(matrix);

  int reduced_nnz = nnz / 2 + SLEQP_MIN(num_cols, num_rows);
  reduced_nnz     = SLEQP_MIN(reduced_nnz, nnz);

  SLEQP_CALL(update_shape(cholmod_data, num_rows, num_cols, reduced_nnz));

  assert(cholmod_data->rhs);
  assert(cholmod_data->rhs->dtype == CHOLMOD_DOUBLE);

  assert(num_cols == num_rows);

  assert(cholmod_data->sparse);
  assert(cholmod_data->sparse->dtype == CHOLMOD_DOUBLE);
  assert(cholmod_data->sparse->itype == CHOLMOD_LONG);

  const double* mat_data = sleqp_sparse_matrix_data(matrix);
  const int* mat_cols    = sleqp_sparse_matrix_cols(matrix);
  const int* mat_rows    = sleqp_sparse_matrix_rows(matrix);

  SuiteSparse_long* cols = cholmod_data->sparse->p;
  SuiteSparse_long* rows = cholmod_data->sparse->i;
  double* data           = cholmod_data->sparse->x;

  int mat_index = 0;
  int col       = 0;

  for (int col = 0; col <= num_cols; ++col)
  {
    cols[col] = 0;
  }

  // extract lower triangular part
  for (int index = 0; index < nnz; ++index)
  {
    while (index >= mat_cols[col + 1])
    {
      ++col;
      cols[col + 1] = cols[col];
    }

    const int row      = mat_rows[index];
    const double value = mat_data[index];

    if (row < col)
    {
      continue;
    }

    data[mat_index] = value;
    rows[mat_index] = row;
    ++(cols[col + 1]);

    ++mat_index;
  }

  assert(mat_index <= reduced_nnz);

  cols[num_cols] = mat_index;

  assert(cholmod_l_check_sparse(cholmod_data->sparse, common) >= 0);

  CHOLMOD_ERROR_CHECK(common);

  cholmod_data->factor = cholmod_l_analyze(cholmod_data->sparse, common);

  CHOLMOD_ERROR_CHECK(common);

  cholmod_l_factorize(cholmod_data->sparse, cholmod_data->factor, common);

  CHOLMOD_ERROR_CHECK(common);

  cholmod_l_print_factor(cholmod_data->factor, "L", common);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
set_cache(double* cache, const SleqpVec* vec)
{
  for (int k = 0; k < vec->nnz; ++k)
  {
    cache[vec->indices[k]] = vec->data[k];
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
reset_cache(double* cache, const SleqpVec* vec)
{
  for (int k = 0; k < vec->nnz; ++k)
  {
    cache[vec->indices[k]] = 0.;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
cholmod_factorization_solve(void* factorization_data, const SleqpVec* rhs)
{
  CHOLMODData* cholmod_data = (CHOLMODData*)factorization_data;

  cholmod_common* common = &(cholmod_data->common);

  SLEQP_CALL(set_cache(cholmod_data->rhs->x, rhs));

  assert(cholmod_l_check_dense(cholmod_data->rhs, common) >= 0);

  if (cholmod_data->sol)
  {
    cholmod_l_free_dense(&cholmod_data->sol, common);
  }

  cholmod_data->sol = cholmod_l_solve(CHOLMOD_A,
                                      cholmod_data->factor,
                                      cholmod_data->rhs,
                                      common);

  assert(cholmod_data->sol);
  assert(cholmod_data->sol->dtype == CHOLMOD_DOUBLE);

  SLEQP_CALL(reset_cache(cholmod_data->rhs->x, rhs));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
cholmod_factorization_condition_estimate(void* factorization_data,
                                         double* condition_estimate)
{
  CHOLMODData* cholmod_data = (CHOLMODData*)factorization_data;

  cholmod_common* common = &(cholmod_data->common);

  const double rcond = cholmod_l_rcond(cholmod_data->factor, common);

  (*condition_estimate) = 1. / rcond;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
cholmod_factorization_solution(void* factorization_data,
                               SleqpVec* sol,
                               int begin,
                               int end,
                               double zero_eps)
{
  CHOLMODData* cholmod_data = (CHOLMODData*)factorization_data;

  assert(begin <= end);

  double* sol_ptr = (double*)cholmod_data->sol->x;

  SLEQP_CALL(
    sleqp_vec_set_from_raw(sol, sol_ptr + begin, end - begin, zero_eps));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
cholmod_factorization_free(void** star)
{
  CHOLMODData* cholmod_data = (CHOLMODData*)(*star);

  if (!cholmod_data)
  {
    return SLEQP_OKAY;
  }

  cholmod_common* common = &(cholmod_data->common);

  cholmod_l_free_dense(&cholmod_data->sol, common);

  cholmod_l_free_dense(&cholmod_data->rhs, common);

  cholmod_l_free_factor(&cholmod_data->factor, common);

  cholmod_l_free_sparse(&cholmod_data->sparse, common);

  cholmod_l_finish(&cholmod_data->common);

  sleqp_free(&cholmod_data);

  *star = NULL;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
cholmod_data_create(CHOLMODData** star)
{
  SLEQP_CALL(sleqp_malloc(star));

  CHOLMODData* cholmod_data = *star;

  *cholmod_data = (CHOLMODData){0};

  cholmod_l_start(&cholmod_data->common);

  cholmod_data->common.error_handler = report_error;
  cholmod_data->common.dtype         = CHOLMOD_DOUBLE;

  cholmod_data->num_cols = SLEQP_NONE;
  cholmod_data->num_rows = SLEQP_NONE;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_factorization_cholmod_create(SleqpFactorization** star,
                                   SleqpParams* params)
{

  SleqpFactorizationCallbacks callbacks
    = {.set_matrix         = cholmod_factorization_set_matrix,
       .solve              = cholmod_factorization_solve,
       .solution           = cholmod_factorization_solution,
       .condition_estimate = cholmod_factorization_condition_estimate,
       .free               = cholmod_factorization_free};

  CHOLMODData* cholmod_data;

  SLEQP_CALL(cholmod_data_create(&cholmod_data));

  SLEQP_CALL(sleqp_factorization_create(star,
                                        SLEQP_FACT_CHOLMOD_NAME,
                                        SLEQP_FACT_CHOLMOD_VERSION,
                                        params,
                                        &callbacks,
                                        CHOLMOD_FLAGS,
                                        (void*)cholmod_data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_factorization_create_default(SleqpFactorization** star,
                                   SleqpParams* params)
{
  SLEQP_CALL(sleqp_factorization_cholmod_create(star, params));

  return SLEQP_OKAY;
}
