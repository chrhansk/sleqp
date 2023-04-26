#include "fact_lapack.h"

#include "mem.h"
#include <assert.h>
#include <string.h>

void
dgetrf_(int* M, int* N, double* A, int* LDA, int* IPIV, int* INFO);

void
dgetrs_(char* TRANS,
        int* N,
        int* NRHS,
        double* A,
        int* LDA,
        int* IPIV,
        double* B,
        int* LDB,
        int* INFO);

typedef struct
{
  SleqpParams* params;

  int max_size;

  int rows;
  int max_rows;

  double* values;
  int* ipiv;

  double* sol;

} LAPACKData;

static SLEQP_RETCODE
lapack_data_create(LAPACKData** star, SleqpParams* params)
{
  SLEQP_CALL(sleqp_malloc(star));

  LAPACKData* lapack_data = *star;

  *lapack_data = (LAPACKData){0};

  SLEQP_CALL(sleqp_params_capture(params));
  lapack_data->params = params;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
store_matrix_values(SleqpMat* matrix, double* values)
{
  const int num_cols = sleqp_mat_num_cols(matrix);

  const int* rows    = sleqp_mat_rows(matrix);
  const int* cols    = sleqp_mat_cols(matrix);
  const double* data = sleqp_mat_data(matrix);

  for (int col = 0; col < num_cols; ++col)
  {
    for (int index = cols[col]; index < cols[col + 1]; ++index)
    {
      const int row = rows[index];

      values[row * num_cols + col] = data[index];
      values[col * num_cols + row] = data[index];
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
lapack_set_matrix(void* fact_data, SleqpMat* matrix)
{
  LAPACKData* lapack_data = (LAPACKData*)fact_data;

  assert(sleqp_mat_num_cols(matrix) == sleqp_mat_num_rows(matrix));

  const int num_rows    = sleqp_mat_num_rows(matrix);
  const int matrix_size = num_rows * num_rows;

  if (lapack_data->max_size < matrix_size)
  {
    SLEQP_CALL(sleqp_realloc(&lapack_data->values, matrix_size));
    lapack_data->max_size = matrix_size;
  }

  if (lapack_data->max_rows < num_rows)
  {
    SLEQP_CALL(sleqp_realloc(&lapack_data->ipiv, num_rows));
    SLEQP_CALL(sleqp_realloc(&lapack_data->sol, num_rows));

    lapack_data->max_rows = num_rows;
  }

  for (int i = 0; i < matrix_size; ++i)
  {
    lapack_data->values[i] = 0.;
  }

  SLEQP_CALL(store_matrix_values(matrix, lapack_data->values));

  lapack_data->rows = num_rows;

  int INFO;

  dgetrf_(&lapack_data->rows,
          &lapack_data->rows,
          lapack_data->values,
          &lapack_data->rows,
          lapack_data->ipiv,
          &INFO);

  if (INFO != 0)
  {
    sleqp_raise(SLEQP_INTERNAL_ERROR, "Failed to factorize using LAPACK");
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
lapack_solve(void* fact_data, const SleqpVec* rhs)
{
  LAPACKData* lapack_data = (LAPACKData*)fact_data;

  assert(rhs->dim == lapack_data->rows);

  SLEQP_CALL(sleqp_vec_to_raw(rhs, lapack_data->sol));

  int one = 1;
  int INFO;
  char* trans = "N";

  dgetrs_(trans,
          &lapack_data->rows,
          &one,
          lapack_data->values,
          &lapack_data->rows,
          lapack_data->ipiv,
          lapack_data->sol,
          &lapack_data->rows,
          &INFO);

  if (INFO != 0)
  {
    sleqp_raise(SLEQP_INTERNAL_ERROR, "Failed to solve using LAPACK");
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
lapack_solution(void* fact_data,
                SleqpVec* sol,
                int begin,
                int end,
                double zero_eps)
{
  LAPACKData* lapack_data = (LAPACKData*)fact_data;

  SLEQP_CALL(sleqp_vec_set_from_raw(sol,
                                    lapack_data->sol + begin,
                                    end - begin,
                                    zero_eps));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
lapack_free(void** star)
{
  LAPACKData* lapack_data = (LAPACKData*)(*star);

  sleqp_free(&lapack_data->sol);
  sleqp_free(&lapack_data->ipiv);
  sleqp_free(&lapack_data->values);

  SLEQP_CALL(sleqp_params_release(&lapack_data->params));

  sleqp_free(&lapack_data);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_fact_lapack_create(SleqpFact** star, SleqpParams* params)
{
  SleqpFactCallbacks callbacks = {.set_matrix = lapack_set_matrix,
                                  .solve      = lapack_solve,
                                  .solution   = lapack_solution,
                                  .condition  = NULL,
                                  .free       = lapack_free};

  LAPACKData* lapack_data = NULL;

  SLEQP_CALL(lapack_data_create(&lapack_data, params));

  SLEQP_CALL(sleqp_fact_create(star,
                               SLEQP_FACT_LAPACK_NAME,
                               SLEQP_FACT_LAPACK_VERSION,
                               params,
                               &callbacks,
                               SLEQP_FACT_FLAGS_LOWER,
                               (void*)lapack_data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_fact_create_default(SleqpFact** star, SleqpParams* params)
{
  SLEQP_CALL(sleqp_fact_lapack_create(star, params));

  return SLEQP_OKAY;
}
