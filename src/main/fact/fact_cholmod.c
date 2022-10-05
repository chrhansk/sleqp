#include "fact_cholmod.h"

#include <cholmod.h>
#include <cholmod_core.h>

#include "cholmod_helpers.h"
#include "defs.h"
#include "error.h"
#include "fail.h"
#include "log.h"
#include "mem.h"

#include "fact/fact.h"

#define CHOLMOD_FLAGS (SLEQP_FACT_FLAGS_PSD | SLEQP_FACT_FLAGS_LOWER)

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

static SLEQP_RETCODE
update_shape(CHOLMODData* cholmod_data, int num_rows, int num_cols, int nnz_max)
{
  cholmod_common* common = &(cholmod_data->common);

  if (num_rows != cholmod_data->num_rows)
  {
    cholmod_l_free_dense(&cholmod_data->rhs, common);

    cholmod_data->rhs
      = cholmod_l_allocate_dense(num_rows, 1, num_rows, CHOLMOD_REAL, common);

    SLEQP_CHOLMOD_ERROR_CHECK(common);

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

    SLEQP_CHOLMOD_ERROR_CHECK(common);
  }
  else
  {
    if (nnz_max > cholmod_data->sparse->nzmax)
    {
      cholmod_l_reallocate_sparse(nnz_max, cholmod_data->sparse, common);

      SLEQP_CHOLMOD_ERROR_CHECK(common);
    }
  }

  cholmod_data->num_cols = num_cols;
  cholmod_data->num_rows = num_rows;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
cholmod_fact_set_matrix(void* fact_data, SleqpSparseMatrix* matrix)
{
  CHOLMODData* cholmod_data = (CHOLMODData*)fact_data;

  cholmod_common* common = &(cholmod_data->common);

  const int num_rows = sleqp_sparse_matrix_num_rows(matrix);
  const int num_cols = sleqp_sparse_matrix_num_cols(matrix);
  const int nnz      = sleqp_sparse_matrix_nnz(matrix);

  SLEQP_CALL(update_shape(cholmod_data, num_rows, num_cols, nnz));

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

  for (int col = 0; col <= num_cols; ++col)
  {
    cols[col] = mat_cols[col];
  }

  // extract lower triangular part
  for (int index = 0; index < nnz; ++index)
  {
    data[index] = mat_data[index];
    rows[index] = mat_rows[index];
  }

  assert(cholmod_l_check_sparse(cholmod_data->sparse, common) >= 0);

  SLEQP_CHOLMOD_ERROR_CHECK(common);

  cholmod_data->factor = cholmod_l_analyze(cholmod_data->sparse, common);

  SLEQP_CHOLMOD_ERROR_CHECK(common);

  cholmod_l_factorize(cholmod_data->sparse, cholmod_data->factor, common);

  SLEQP_CHOLMOD_ERROR_CHECK(common);

  // cholmod_l_print_factor(cholmod_data->factor, "L", common);

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
cholmod_fact_solve(void* fact_data, const SleqpVec* rhs)
{
  CHOLMODData* cholmod_data = (CHOLMODData*)fact_data;

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
cholmod_fact_condition(void* fact_data, double* condition)
{
  CHOLMODData* cholmod_data = (CHOLMODData*)fact_data;

  cholmod_common* common = &(cholmod_data->common);

  const double rcond = cholmod_l_rcond(cholmod_data->factor, common);

  (*condition) = 1. / rcond;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
cholmod_fact_solution(void* fact_data,
                      SleqpVec* sol,
                      int begin,
                      int end,
                      double zero_eps)
{
  CHOLMODData* cholmod_data = (CHOLMODData*)fact_data;

  assert(begin <= end);

  double* sol_ptr = (double*)cholmod_data->sol->x;

  SLEQP_CALL(
    sleqp_vec_set_from_raw(sol, sol_ptr + begin, end - begin, zero_eps));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
cholmod_fact_free(void** star)
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

  cholmod_data->common.error_handler = sleqp_cholmod_report_error;
  cholmod_data->common.dtype         = CHOLMOD_DOUBLE;

  cholmod_data->num_cols = SLEQP_NONE;
  cholmod_data->num_rows = SLEQP_NONE;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_fact_cholmod_create(SleqpFact** star, SleqpParams* params)
{

  SleqpFactorizationCallbacks callbacks
    = {.set_matrix = cholmod_fact_set_matrix,
       .solve      = cholmod_fact_solve,
       .solution   = cholmod_fact_solution,
       .condition  = cholmod_fact_condition,
       .free       = cholmod_fact_free};

  CHOLMODData* cholmod_data;

  SLEQP_CALL(cholmod_data_create(&cholmod_data));

  SLEQP_CALL(sleqp_fact_create(star,
                               SLEQP_FACT_CHOLMOD_NAME,
                               SLEQP_FACT_CHOLMOD_VERSION,
                               params,
                               &callbacks,
                               CHOLMOD_FLAGS,
                               (void*)cholmod_data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_fact_create_default(SleqpFact** star, SleqpParams* params)
{
  SLEQP_CALL(sleqp_fact_cholmod_create(star, params));

  return SLEQP_OKAY;
}
