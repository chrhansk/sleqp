#include "hsl_matrix.h"

#include <assert.h>

#include "sleqp_mem.h"

SLEQP_RETCODE hsl_matrix_set(HSLMatrix* hsl_matrix,
                             SleqpSparseMatrix* matrix)
{
  const int num_rows = sleqp_sparse_matrix_get_num_rows(matrix);
  const int num_cols = sleqp_sparse_matrix_get_num_cols(matrix);

  assert(num_rows == num_cols);

  hsl_matrix->dim = num_rows;

  const double* matrix_data = sleqp_sparse_matrix_get_data(matrix);
  const int* matrix_cols = sleqp_sparse_matrix_get_cols(matrix);
  const int* matrix_rows = sleqp_sparse_matrix_get_rows(matrix);
  const int matrix_nnz = sleqp_sparse_matrix_get_nnz(matrix);

  if(hsl_matrix->max_nnz < matrix_nnz)
  {
    SLEQP_CALL(sleqp_realloc(&(hsl_matrix->cols), matrix_nnz));
    SLEQP_CALL(sleqp_realloc(&(hsl_matrix->rows), matrix_nnz));
    SLEQP_CALL(sleqp_realloc(&(hsl_matrix->data), matrix_nnz));

    hsl_matrix->max_nnz = matrix_nnz;
  }

  double* hsl_data = hsl_matrix->data;
  int32_t* hsl_rows = hsl_matrix->rows;
  int32_t* hsl_cols = hsl_matrix->cols;

  int32_t hsl_pos = 0;
  int32_t col = 0;

  for(int index = 0; index < matrix_nnz; ++index)
  {
    while(index >= matrix_cols[col + 1])
    {
      ++col;
    }

    const int32_t row = matrix_rows[index];
    const double entry = matrix_data[index];

    // Convert indices to be 1-based
    if(row <= col)
    {
      hsl_data[hsl_pos] = entry;
      hsl_cols[hsl_pos] = col + 1;
      hsl_rows[hsl_pos] = row + 1;

      ++hsl_pos;
    }

    // TODO: Skip the rest of the column
  }

  hsl_matrix->nnz = hsl_pos;

  return SLEQP_OKAY;
}

SLEQP_RETCODE hsl_matrix_clear(HSLMatrix* hsl_matrix)
{
  sleqp_free(&(hsl_matrix->cols));
  sleqp_free(&(hsl_matrix->rows));
  sleqp_free(&(hsl_matrix->data));

  hsl_matrix->nnz = 0;
  hsl_matrix->max_nnz = 0;

  return SLEQP_OKAY;
}
