#include "sleqp_sparse_factorization_mumps.h"

#include <assert.h>

#include <dmumps_c.h>

#include "sleqp_mem.h"

static const bool sleqp_mumps_verbose = true;

typedef struct SleqpMUMPSData
{
	DMUMPS_STRUC_C id;

  int dim;
  int nnz;
  int max_nnz;

  int* cols;
  int* rows;
  double* data;

  double* rhs_sol;

} SleqpMUMPSData;

#define SLEQP_MUMPS_IS_ERROR(value) (value < 0)

#define SLEQP_MUMPS_CALL(id)                                            \
  do                                                                    \
  {                                                                     \
    dmumps_c(&(id));                                                    \
    int sleqp_mumps_status = ((id).infog[0]);                           \
                                                                        \
    if(SLEQP_MUMPS_IS_ERROR(sleqp_mumps_status))                        \
    {                                                                   \
                                                                        \
      sleqp_log_error("Caught MUMPS error <%d> in function %s",         \
                      sleqp_mumps_status,                               \
                      __func__);                                        \
                                                                        \
      return SLEQP_INTERNAL_ERROR;                                      \
    }                                                                   \
  } while(0)


static SLEQP_RETCODE sleqp_mumps_create(SleqpMUMPSData** star)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpMUMPSData* sleqp_mumps_data = *star;

  *sleqp_mumps_data = (SleqpMUMPSData){0};

  if(sleqp_mumps_verbose)
  {
    sleqp_mumps_data->id.icntl[3] = 3;
  }
  else
  {
    sleqp_mumps_data->id.icntl[3] = 1;
  }

  // Mark as symmetric
  sleqp_mumps_data->id.sym = 2;

  // Use local parallelization
  sleqp_mumps_data->id.par = 1;

  // init job
  sleqp_mumps_data->id.job = -1;
  SLEQP_MUMPS_CALL(sleqp_mumps_data->id);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE matrix_reserve(SleqpMUMPSData* sleqp_mumps_data,
                                    int nnz)
{
  if(sleqp_mumps_data->max_nnz < nnz)
  {
    SLEQP_CALL(sleqp_realloc(&(sleqp_mumps_data->cols), nnz));
    SLEQP_CALL(sleqp_realloc(&(sleqp_mumps_data->rows), nnz));
    SLEQP_CALL(sleqp_realloc(&(sleqp_mumps_data->data), nnz));

    sleqp_mumps_data->max_nnz = nnz;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE matrix_fill(SleqpMUMPSData* sleqp_mumps_data,
                                 SleqpSparseMatrix* matrix)
{
  const int nnz = sleqp_sparse_matrix_get_nnz(matrix);

  const int num_rows = sleqp_sparse_matrix_get_num_rows(matrix);
  const int num_cols = sleqp_sparse_matrix_get_num_cols(matrix);

  assert(num_rows == num_cols);

  const int dim = num_rows;

  if(sleqp_mumps_data->dim < dim)
  {
    SLEQP_CALL(sleqp_realloc(&(sleqp_mumps_data->rhs_sol), dim));
  }

  SLEQP_CALL(matrix_reserve(sleqp_mumps_data, nnz));


  sleqp_mumps_data->dim = dim;

  const double* data = sleqp_sparse_matrix_get_data(matrix);
  const int* cols = sleqp_sparse_matrix_get_cols(matrix);
  const int* rows = sleqp_sparse_matrix_get_rows(matrix);

  int col = 0;

  int data_index = 0;

  for(int index = 0; index < nnz; ++index)
  {
    while(index >= cols[col + 1])
    {
      ++col;
    }

    const int row = rows[index];
    const double value = data[index];

    if(row < col)
    {
      continue;
    }

    sleqp_mumps_data->rows[data_index] = (row + 1);
    sleqp_mumps_data->cols[data_index] = (col + 1);
    sleqp_mumps_data->data[data_index] = value;

    ++data_index;
  }

  assert(data_index <= nnz);

  sleqp_mumps_data->nnz = data_index;

  return SLEQP_OKAY;
}


static SLEQP_RETCODE sleqp_mumps_set_matrix(void* factorization_data,
                                            SleqpSparseMatrix* matrix)
{
  SleqpMUMPSData* sleqp_mumps_data = (SleqpMUMPSData*) factorization_data;

  SLEQP_CALL(matrix_fill(sleqp_mumps_data, matrix));

  sleqp_mumps_data->id.n = sleqp_mumps_data->dim;
	sleqp_mumps_data->id.nz = sleqp_mumps_data->nnz;

  sleqp_mumps_data->id.irn = sleqp_mumps_data->rows;
	sleqp_mumps_data->id.jcn = sleqp_mumps_data->cols;
	sleqp_mumps_data->id.a = sleqp_mumps_data->data;

  // analysis job
  sleqp_mumps_data->id.job = 1;
  SLEQP_MUMPS_CALL(sleqp_mumps_data->id);

  // factorization job
  sleqp_mumps_data->id.job = 2;
  SLEQP_MUMPS_CALL(sleqp_mumps_data->id);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE sleqp_mumps_solve(void* factorization_data,
                                       SleqpSparseVec* rhs)
{
  SleqpMUMPSData* sleqp_mumps_data = (SleqpMUMPSData*) factorization_data;

  const int dim = sleqp_mumps_data->dim;

  assert(rhs->dim == dim);

  SLEQP_CALL(sleqp_sparse_vector_to_raw(rhs, sleqp_mumps_data->rhs_sol));

	sleqp_mumps_data->id.rhs = sleqp_mumps_data->rhs_sol;
	sleqp_mumps_data->id.nrhs = 1;
	sleqp_mumps_data->id.lrhs = dim;

  // substitution job
  sleqp_mumps_data->id.job = 3;
  SLEQP_MUMPS_CALL(sleqp_mumps_data->id);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE sleqp_mumps_get_sol(void* factorization_data,
                                         SleqpSparseVec* sol,
                                         int begin,
                                         int end,
                                         double zero_eps)
{
  SleqpMUMPSData* sleqp_mumps_data = (SleqpMUMPSData*) factorization_data;

  SLEQP_CALL(sleqp_sparse_vector_from_raw(sol,
                                          sleqp_mumps_data->rhs_sol + begin,
                                          end - begin,
                                          zero_eps));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE sleqp_mumps_get_condition_estimate(void* factorization_data,
                                                        double* condition_estimate)
{
  //SleqpMUMPSData* sleqp_mumps_data = (SleqpMUMPSData*) factorization_data;

  (*condition_estimate) = SLEQP_NONE;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE sleqp_mumps_free(void **star)
{
  SleqpMUMPSData* sleqp_mumps_data = (SleqpMUMPSData*) (*star);

  // de-init job
  sleqp_mumps_data->id.job = -2;
  SLEQP_MUMPS_CALL(sleqp_mumps_data->id);

  sleqp_free(&sleqp_mumps_data->cols);
  sleqp_free(&sleqp_mumps_data->rows);
  sleqp_free(&sleqp_mumps_data->data);
  sleqp_free(&sleqp_mumps_data->rhs_sol);

  sleqp_free(&sleqp_mumps_data);

  *star = NULL;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_factorization_mumps_create(SleqpSparseFactorization** star,
                                                      SleqpParams* params)
{
  SleqpSparseFactorizationCallbacks callbacks = {
    .set_matrix = sleqp_mumps_set_matrix,
    .solve = sleqp_mumps_solve,
    .get_sol = sleqp_mumps_get_sol,
    .get_condition_estimate = sleqp_mumps_get_condition_estimate,
    .free = sleqp_mumps_free
  };

  SleqpMUMPSData* sleqp_mumps_data;

  SLEQP_CALL(sleqp_mumps_create(&sleqp_mumps_data));

  SLEQP_CALL(sleqp_sparse_factorization_create(star,
                                               params,
                                               &callbacks,
                                               (void*) sleqp_mumps_data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_factorization_create_default(SleqpSparseFactorization** star,
                                                        SleqpParams* params)
{
  SLEQP_CALL(sleqp_sparse_factorization_mumps_create(star, params));

  return SLEQP_OKAY;
}
