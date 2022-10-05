#include "fact_spqr.h"

#include <SuiteSparseQR_C.h>
#include <cholmod.h>

#include "cholmod_helpers.h"
#include "defs.h"
#include "error.h"
#include "fail.h"
#include "log.h"
#include "mem.h"

typedef struct SPQRData
{
  SleqpSparseMatrix* matrix;

  cholmod_common common;
  cholmod_sparse* sparse;

  SuiteSparseQR_C_factorization* factorization;

  cholmod_dense* sol;
  cholmod_dense* rhs;

  int size;

} SPQRData;

static SLEQP_RETCODE
update_shape(SPQRData* spqr, int size, int nnz_max)
{
  cholmod_common* common = &(spqr->common);

  if (size != spqr->size)
  {
    cholmod_l_free_dense(&spqr->rhs, common);

    spqr->rhs = cholmod_l_allocate_dense(size, 1, size, CHOLMOD_REAL, common);

    SLEQP_CHOLMOD_ERROR_CHECK(common);

    for (int i = 0; i < size; ++i)
    {
      ((double*)spqr->rhs->x)[i] = 0.;
    }
  }

  if ((!spqr->sparse) || (size != spqr->size))
  {
    cholmod_l_free_sparse(&spqr->sparse, common);

    spqr->sparse
      = cholmod_l_allocate_sparse(size,
                                  size,
                                  nnz_max,
                                  true, // sorted
                                  true, // packed
                                  0,    // both upper/lower parts are stored
                                  CHOLMOD_REAL,
                                  common);

    SLEQP_CHOLMOD_ERROR_CHECK(common);
  }
  else
  {
    if (nnz_max > spqr->sparse->nzmax)
    {
      cholmod_l_reallocate_sparse(nnz_max, spqr->sparse, common);

      SLEQP_CHOLMOD_ERROR_CHECK(common);
    }
  }

  spqr->size = size;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
spqr_factorization_set_matrix(void* fact_data, SleqpSparseMatrix* matrix)
{
  SPQRData* spqr = (SPQRData*)fact_data;

  cholmod_common* common = &(spqr->common);

  const int num_rows = sleqp_sparse_matrix_num_rows(matrix);
  const int num_cols = sleqp_sparse_matrix_num_cols(matrix);
  const int nnz      = sleqp_sparse_matrix_nnz(matrix);

  assert(num_rows == num_cols);

  SLEQP_CALL(update_shape(spqr, num_rows, nnz));

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

  // cholmod_l_write_sparse(stdout, spqr->sparse, NULL, NULL, common);

  spqr->factorization = SuiteSparseQR_C_factorize(SPQR_ORDERING_DEFAULT,
                                                  SPQR_DEFAULT_TOL,
                                                  spqr->sparse,
                                                  common);

  SLEQP_CHOLMOD_ERROR_CHECK(common);

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
spqr_factorization_solve(void* fact_data, const SleqpVec* rhs)
{
  SPQRData* spqr = (SPQRData*)fact_data;

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

  SLEQP_CHOLMOD_ERROR_CHECK(common);

  // X = R\(E*Y)
  spqr->sol
    = SuiteSparseQR_C_solve(SPQR_RETX_EQUALS_B, spqr->factorization, y, common);

  SLEQP_CHOLMOD_ERROR_CHECK(common);

  cholmod_l_free_dense(&y, common);

  assert(spqr->sol);
  assert(spqr->sol->dtype == CHOLMOD_DOUBLE);

  SLEQP_CALL(reset_cache(spqr->rhs->x, rhs));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
spqr_factorization_solution(void* fact_data,
                            SleqpVec* sol,
                            int begin,
                            int end,
                            double zero_eps)
{
  SPQRData* spqr = (SPQRData*)fact_data;

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

  spqr->common.error_handler = sleqp_cholmod_report_error;
  spqr->common.dtype         = CHOLMOD_DOUBLE;
  spqr->size                 = SLEQP_NONE;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_fact_spqr_create(SleqpFact** star, SleqpParams* params)
{

  SleqpFactorizationCallbacks callbacks
    = {.set_matrix = spqr_factorization_set_matrix,
       .solve      = spqr_factorization_solve,
       .solution   = spqr_factorization_solution,
       .condition  = NULL,
       .free       = spqr_factorization_free};

  SPQRData* spqr_data;

  SLEQP_CALL(spqr_data_create(&spqr_data));

  SLEQP_CALL(sleqp_fact_create(star,
                               SLEQP_FACT_SPQR_NAME,
                               SLEQP_FACT_SPQR_VERSION,
                               params,
                               &callbacks,
                               SLEQP_FACT_FLAGS_NONE,
                               (void*)spqr_data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_fact_create_default(SleqpFact** star, SleqpParams* params)
{
  SLEQP_CALL(sleqp_fact_spqr_create(star, params));

  return SLEQP_OKAY;
}
