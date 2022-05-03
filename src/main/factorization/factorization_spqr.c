#include "factorization_spqr.h"

#include "fail.h"
#include <assert.h>

#include <SuiteSparseQR_C.h>
#include <cholmod.h>

#include "defs.h"
#include "error.h"
#include "log.h"
#include "mem.h"
#include "pub_log.h"
#include "pub_types.h"
#include "sparse/pub_sparse_matrix.h"

typedef struct SPQRData
{
  SleqpSparseMatrix* matrix;

  cholmod_common common;
  cholmod_sparse* sparse;

  SuiteSparseQR_C_factorization* factorization;

  cholmod_dense* sol;
  cholmod_dense* rhs;

  int num_rows;
  int num_cols;

} SPQRData;

static void
report_error(int status, const char* file, int line, const char* message)
{
  sleqp_log_error("CHOLMOD error %s", message);
}

static SLEQP_RETCODE
spqr_error_string(int status, const char** message)
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

#define SPQR_ERROR_CHECK(common)                                               \
  do                                                                           \
  {                                                                            \
    if ((common)->status < 0)                                                  \
    {                                                                          \
      const char* message;                                                     \
      SLEQP_CALL(spqr_error_string((common)->status, &message));               \
      sleqp_raise(SLEQP_INTERNAL_ERROR,                                        \
                  "Encountered error in SPQR: %s",                             \
                  message);                                                    \
    }                                                                          \
  } while (false)

static SLEQP_RETCODE
update_shape(SPQRData* spqr, int num_rows, int num_cols, int nnz_max)
{
  cholmod_common* common = &(spqr->common);

  if (num_rows != spqr->num_rows)
  {
    cholmod_l_free_dense(&spqr->rhs, common);

    spqr->rhs
      = cholmod_l_allocate_dense(num_rows, 1, num_rows, CHOLMOD_REAL, common);

    SPQR_ERROR_CHECK(common);

    for (int i = 0; i < num_rows; ++i)
    {
      ((double*)spqr->rhs->x)[i] = 0.;
    }
  }

  if ((!spqr->sparse) || (num_cols != spqr->num_cols)
      || (num_rows != spqr->num_rows))
  {
    cholmod_l_free_sparse(&spqr->sparse, common);

    spqr->sparse
      = cholmod_l_allocate_sparse(num_rows,
                                  num_cols,
                                  nnz_max,
                                  true, // sorted
                                  true, // packed
                                  0,    // both upper/lower parts are stored
                                  CHOLMOD_REAL,
                                  common);

    SPQR_ERROR_CHECK(common);
  }
  else
  {
    if (nnz_max > spqr->sparse->nzmax)
    {
      cholmod_l_reallocate_sparse(nnz_max, spqr->sparse, common);

      SPQR_ERROR_CHECK(common);
    }
  }

  spqr->num_cols = num_cols;
  spqr->num_rows = num_rows;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
spqr_factorization_set_matrix(void* factorization_data,
                              SleqpSparseMatrix* matrix)
{
  SPQRData* spqr = (SPQRData*)factorization_data;

  cholmod_common* common = &(spqr->common);

  const int num_rows = sleqp_sparse_matrix_num_rows(matrix);
  const int num_cols = sleqp_sparse_matrix_num_cols(matrix);
  const int nnz      = sleqp_sparse_matrix_nnz(matrix);

  SLEQP_CALL(update_shape(spqr, num_rows, num_cols, nnz));

  assert(spqr->rhs);
  assert(spqr->rhs->dtype == CHOLMOD_DOUBLE);

  assert(num_cols == num_rows);

  assert(spqr->sparse);
  assert(spqr->sparse->dtype == CHOLMOD_DOUBLE);
  assert(spqr->sparse->itype == CHOLMOD_LONG);

  const double* mat_data = sleqp_sparse_matrix_data(matrix);
  const int* mat_cols    = sleqp_sparse_matrix_cols(matrix);
  const int* mat_rows    = sleqp_sparse_matrix_rows(matrix);

  SuiteSparse_long* cols = spqr->sparse->p;
  SuiteSparse_long* rows = spqr->sparse->i;
  double* data           = spqr->sparse->x;

  for (int i = 0; i < nnz; ++i)
  {
    data[i] = mat_data[i];
  }

  for (int i = 0; i <= num_cols; ++i)
  {
    cols[i] = mat_cols[i];
  }

  for (int i = 0; i < nnz; ++i)
  {
    rows[i] = mat_rows[i];
  }

  assert(cholmod_l_check_sparse(spqr->sparse, common) >= 0);

  spqr->factorization = SuiteSparseQR_C_factorize(SPQR_ORDERING_DEFAULT,
                                                  SPQR_DEFAULT_TOL,
                                                  spqr->sparse,
                                                  common);

  SPQR_ERROR_CHECK(common);

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
spqr_factorization_solve(void* factorization_data, const SleqpVec* rhs)
{
  SPQRData* spqr = (SPQRData*)factorization_data;

  cholmod_common* common = &(spqr->common);

  SLEQP_CALL(set_cache(spqr->rhs->x, rhs));

  assert(cholmod_l_check_dense(spqr->rhs, common) >= 0);

  if (spqr->sol)
  {
    cholmod_l_free_dense(&spqr->sol, common);
  }

  // Y = Q'*B
  cholmod_dense* y
    = SuiteSparseQR_C_qmult(SPQR_QTX, spqr->factorization, spqr->rhs, common);

  SPQR_ERROR_CHECK(common);

  // X = R\(E*Y)
  spqr->sol
    = SuiteSparseQR_C_solve(SPQR_RETX_EQUALS_B, spqr->factorization, y, common);

  SPQR_ERROR_CHECK(common);

  cholmod_l_free_dense(&y, common);

  assert(spqr->sol);
  assert(spqr->sol->dtype == CHOLMOD_DOUBLE);

  SLEQP_CALL(reset_cache(spqr->rhs->x, rhs));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
spqr_factorization_condition_estimate(void* factorization_data,
                                      double* condition_estimate)
{
  (*condition_estimate) = SLEQP_NONE;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
spqr_factorization_solution(void* factorization_data,
                            SleqpVec* sol,
                            int begin,
                            int end,
                            double zero_eps)
{
  SPQRData* spqr = (SPQRData*)factorization_data;

  assert(begin <= end);

  double* sol_ptr = (double*)spqr->sol->x;

  SLEQP_CALL(
    sleqp_vec_set_from_raw(sol, sol_ptr + begin, end - begin, zero_eps));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
spqr_factorization_free(void** star)
{
  SPQRData* spqr = (SPQRData*)(*star);

  if (!spqr)
  {
    return SLEQP_OKAY;
  }

  cholmod_common* common = &(spqr->common);

  cholmod_l_free_dense(&spqr->sol, common);

  cholmod_l_free_dense(&spqr->rhs, common);

  SuiteSparseQR_C_free(&spqr->factorization, common);

  cholmod_l_free_sparse(&spqr->sparse, common);

  cholmod_l_finish(&spqr->common);

  sleqp_free(&spqr);

  *star = NULL;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
spqr_data_create(SPQRData** star)
{
  SLEQP_CALL(sleqp_malloc(star));

  SPQRData* spqr = *star;

  *spqr = (SPQRData){0};

  cholmod_l_start(&spqr->common);

  spqr->common.error_handler = report_error;
  spqr->common.dtype         = CHOLMOD_DOUBLE;

  spqr->num_cols = SLEQP_NONE;
  spqr->num_rows = SLEQP_NONE;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_factorization_spqr_create(SleqpFactorization** star, SleqpParams* params)
{

  SleqpFactorizationCallbacks callbacks
    = {.set_matrix         = spqr_factorization_set_matrix,
       .solve              = spqr_factorization_solve,
       .solution           = spqr_factorization_solution,
       .condition_estimate = spqr_factorization_condition_estimate,
       .free               = spqr_factorization_free};

  SPQRData* spqr_data;

  SLEQP_CALL(spqr_data_create(&spqr_data));

  SLEQP_CALL(sleqp_factorization_create(star,
                                        SLEQP_FACT_SPQR_NAME,
                                        SLEQP_FACT_SPQR_VERSION,
                                        params,
                                        &callbacks,
                                        (void*)spqr_data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_factorization_create_default(SleqpFactorization** star,
                                   SleqpParams* params)
{
  SLEQP_CALL(sleqp_factorization_spqr_create(star, params));

  return SLEQP_OKAY;
}
