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

  int num_rows;
  int num_cols;

} SPQRData;

static SLEQP_RETCODE
update_shape(SPQRData* spqr, int num_rows, int num_cols, int nnz_max)
{
  cholmod_common* common = &(spqr->common);

  const int rhs_size = SLEQP_MAX(num_rows, num_cols);

  if (!spqr->rhs || (rhs_size != spqr->rhs->nrow))
  {
    cholmod_l_free_dense(&spqr->rhs, common);

    spqr->rhs
      = cholmod_l_allocate_dense(rhs_size, 1, rhs_size, CHOLMOD_REAL, common);

    SLEQP_CHOLMOD_ERROR_CHECK(common);

    for (int i = 0; i < rhs_size; ++i)
    {
      ((double*)spqr->rhs->x)[i] = 0.;
    }
  }

  if ((!spqr->sparse) || (num_rows != spqr->num_rows)
      || (num_cols != spqr->num_cols))
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

  spqr->num_rows = num_rows;
  spqr->num_cols = num_cols;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
spqr_fact_set_matrix(void* fact_data, SleqpSparseMatrix* matrix)
{
  SPQRData* spqr = (SPQRData*)fact_data;

  cholmod_common* common = &(spqr->common);

  const int num_rows = sleqp_sparse_matrix_num_rows(matrix);
  const int num_cols = sleqp_sparse_matrix_num_cols(matrix);
  const int nnz      = sleqp_sparse_matrix_nnz(matrix);

  assert(num_rows >= num_cols);

  SLEQP_CALL(update_shape(spqr, num_rows, num_cols, nnz));

  assert(spqr->rhs);
  assert(spqr->rhs->dtype == CHOLMOD_DOUBLE);

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

  SuiteSparseQR_C_free(&spqr->factorization, common);

  spqr->factorization = SuiteSparseQR_C_factorize(SPQR_ORDERING_DEFAULT,
                                                  SPQR_DEFAULT_TOL,
                                                  spqr->sparse,
                                                  common);

  SLEQP_CHOLMOD_ERROR_CHECK(common);
  SLEQP_CHOLMOD_NULL_CHECK(spqr->factorization);

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
spqr_set_sol(SPQRData* spqr, cholmod_dense* sol)
{
  cholmod_common* common = &(spqr->common);

  SLEQP_CHOLMOD_ERROR_CHECK(common);
  SLEQP_CHOLMOD_NULL_CHECK(sol);

  assert(sol->dtype == CHOLMOD_DOUBLE);

  if (spqr->sol)
  {
    cholmod_l_free_dense(&spqr->sol, common);
  }

  spqr->sol = sol;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
spqr_fact_solve(void* fact_data, const SleqpVec* rhs)
{
  SPQRData* spqr = (SPQRData*)fact_data;

  cholmod_common* common = &(spqr->common);

  spqr->rhs->nrow = rhs->dim;
  SLEQP_CALL(set_cache(spqr->rhs->x, rhs));

  assert(cholmod_l_check_dense(spqr->rhs, common) >= 0);
  // Y = Q'*B
  cholmod_dense* y
    = SuiteSparseQR_C_qmult(SPQR_QTX, spqr->factorization, spqr->rhs, common);

  SLEQP_CHOLMOD_ERROR_CHECK(common);
  SLEQP_CHOLMOD_NULL_CHECK(y);

  // X = R\(E*Y)
  spqr_set_sol(
    spqr,
    SuiteSparseQR_C_solve(SPQR_RETX_EQUALS_B, spqr->factorization, y, common));

  cholmod_l_free_dense(&y, common);

  SLEQP_CALL(reset_cache(spqr->rhs->x, rhs));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
spqr_qr_solve_tri(void* fact_data, const SleqpVec* rhs, SleqpVec* sol)
{
  SPQRData* spqr = (SPQRData*)fact_data;

  cholmod_common* common = &(spqr->common);

  // solve against enlarged rhs with artificial zeros
  // (discarded by SPQR)
  spqr->rhs->nrow = spqr->num_rows;
  SLEQP_CALL(set_cache(spqr->rhs->x, rhs));

  spqr_set_sol(spqr,
               SuiteSparseQR_C_solve(SPQR_RX_EQUALS_B,
                                     spqr->factorization,
                                     spqr->rhs,
                                     common));

  SLEQP_CALL(reset_cache(spqr->rhs->x, rhs));

  double* sol_ptr = (double*)spqr->sol->x;

  SLEQP_CALL(sleqp_vec_set_from_raw(sol, sol_ptr, sol->dim, 0.));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
spqr_qr_solve_tri_trans(void* fact_data, const SleqpVec* rhs, SleqpVec* sol)
{
  SPQRData* spqr = (SPQRData*)fact_data;

  cholmod_common* common = &(spqr->common);

  spqr->rhs->nrow = rhs->dim;
  SLEQP_CALL(set_cache(spqr->rhs->x, rhs));

  spqr_set_sol(spqr,
               SuiteSparseQR_C_solve(SPQR_RTX_EQUALS_B,
                                     spqr->factorization,
                                     spqr->rhs,
                                     common));

  SLEQP_CALL(reset_cache(spqr->rhs->x, rhs));

  double* sol_ptr = (double*)spqr->sol->x;

  SLEQP_CALL(sleqp_vec_set_from_raw(sol, sol_ptr, sol->dim, 0.));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
spqr_qr_mult_orth(void* fact_data, const SleqpVec* direction, SleqpVec* product)
{
  SPQRData* spqr = (SPQRData*)fact_data;

  cholmod_common* common = &(spqr->common);

  spqr->rhs->nrow = spqr->sparse->nrow;
  SLEQP_CALL(set_cache(spqr->rhs->x, direction));

  assert(cholmod_l_check_dense(spqr->rhs, common) >= 0);

  spqr_set_sol(
    spqr,
    SuiteSparseQR_C_qmult(SPQR_QX, spqr->factorization, spqr->rhs, common));

  SLEQP_CALL(reset_cache(spqr->rhs->x, direction));

  double* prod_ptr = (double*)spqr->sol->x;

  SLEQP_CALL(sleqp_vec_set_from_raw(product, prod_ptr, product->dim, 0.));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
spqr_qr_mult_orth_trans(void* fact_data,
                        const SleqpVec* direction,
                        SleqpVec* product)
{
  SPQRData* spqr = (SPQRData*)fact_data;

  cholmod_common* common = &(spqr->common);

  spqr->rhs->nrow = direction->dim;
  SLEQP_CALL(set_cache(spqr->rhs->x, direction));

  assert(cholmod_l_check_dense(spqr->rhs, common) >= 0);

  spqr_set_sol(
    spqr,
    SuiteSparseQR_C_qmult(SPQR_QTX, spqr->factorization, spqr->rhs, common));

  SLEQP_CHOLMOD_NULL_CHECK(spqr->sol);

  SLEQP_CALL(reset_cache(spqr->rhs->x, direction));

  double* prod_ptr = (double*)spqr->sol->x;

  SLEQP_CALL(sleqp_vec_set_from_raw(product, prod_ptr, product->dim, 0.));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
spqr_fact_sol(void* fact_data,
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
spqr_fact_free(void** star)
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
  spqr->num_rows             = SLEQP_NONE;
  spqr->num_cols             = SLEQP_NONE;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_fact_spqr_create(SleqpFactQR** star, SleqpParams* params)
{
  SleqpQRCallbacks callbacks = {.set_matrix      = spqr_fact_set_matrix,
                                .solve_tri       = spqr_qr_solve_tri,
                                .solve_tri_trans = spqr_qr_solve_tri_trans,
                                .mult_orth       = spqr_qr_mult_orth,
                                .mult_orth_trans = spqr_qr_mult_orth_trans,
                                .free            = spqr_fact_free};

  SPQRData* spqr_data;

  SLEQP_CALL(spqr_data_create(&spqr_data));

  SLEQP_CALL(sleqp_qr_create(star,
                             SLEQP_FACT_SPQR_NAME,
                             SLEQP_FACT_SPQR_VERSION,
                             params,
                             &callbacks,
                             (void*)spqr_data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_fact_create_default(SleqpFact** star, SleqpParams* params)
{
  SleqpFactCallbacks callbacks = {.set_matrix = spqr_fact_set_matrix,
                                  .solve      = spqr_fact_solve,
                                  .solution   = spqr_fact_sol,
                                  .free       = spqr_fact_free};

  SPQRData* spqr_data;

  SLEQP_CALL(spqr_data_create(&spqr_data));

  SLEQP_CALL(sleqp_fact_create(star,
                               SLEQP_FACT_SPQR_NAME,
                               SLEQP_FACT_SPQR_VERSION,
                               params,
                               &callbacks,
                               SLEQP_FACT_FLAGS_NONE,
                               spqr_data));

  return SLEQP_OKAY;
}
