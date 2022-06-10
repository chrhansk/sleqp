#include "hsl_matrix.h"

#include <assert.h>

#include "types.h"

#include "log.h"
#include "mem.h"
#include "sparse/sparse_matrix.h"

SLEQP_RETCODE
hsl_matrix_set(HSLMatrix* hsl_matrix, SleqpSparseMatrix* matrix)
{
  const int num_rows = sleqp_sparse_matrix_num_rows(matrix);
  const int num_cols = sleqp_sparse_matrix_num_cols(matrix);

  assert(num_rows == num_cols);

  hsl_matrix->dim = num_rows;

  const double* matrix_data = sleqp_sparse_matrix_data(matrix);
  const int* matrix_cols    = sleqp_sparse_matrix_cols(matrix);
  const int* matrix_rows    = sleqp_sparse_matrix_rows(matrix);
  const int matrix_nnz      = sleqp_sparse_matrix_nnz(matrix);

  if (hsl_matrix->max_nnz < matrix_nnz)
  {
    SLEQP_CALL(sleqp_realloc(&(hsl_matrix->cols), matrix_nnz));
    SLEQP_CALL(sleqp_realloc(&(hsl_matrix->rows), matrix_nnz));
    SLEQP_CALL(sleqp_realloc(&(hsl_matrix->data), matrix_nnz));

    hsl_matrix->max_nnz = matrix_nnz;
  }

  double* hsl_data  = hsl_matrix->data;
  int32_t* hsl_rows = hsl_matrix->rows;
  int32_t* hsl_cols = hsl_matrix->cols;

  for (int index = 0; index < matrix_nnz; ++index)
  {
    hsl_data[index] = matrix_data[index];
  }

  for (int index = 0; index < matrix_nnz; ++index)
  {
    hsl_rows[index] = matrix_rows[index] + 1;
  }

  for (int col = 0; col < num_cols; ++col)
  {
    for (int index = matrix_cols[col]; index < matrix_cols[col + 1]; ++index)
    {
      hsl_cols[index] = col + 1;
    }
  }

  hsl_matrix->nnz = matrix_nnz;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
hsl_matrix_clear(HSLMatrix* hsl_matrix)
{
  sleqp_free(&(hsl_matrix->cols));
  sleqp_free(&(hsl_matrix->rows));
  sleqp_free(&(hsl_matrix->data));

  hsl_matrix->nnz     = 0;
  hsl_matrix->max_nnz = 0;

  return SLEQP_OKAY;
}
