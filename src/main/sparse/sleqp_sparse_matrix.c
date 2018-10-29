#include "sleqp_sparse_matrix.h"

#include <assert.h>

#include "sleqp_cmp.h"
#include "sleqp_mem.h"

SLEQP_RETCODE sleqp_sparse_matrix_create(SleqpSparseMatrix** mstar,
                                         int num_rows,
                                         int num_cols,
                                         int nnz_max)
{
  SLEQP_CALL(sleqp_malloc(mstar));

  SleqpSparseMatrix* matrix = *mstar;

  matrix->nnz = 0;
  matrix->nnz_max = nnz_max;

  matrix->num_cols = num_cols;
  matrix->num_rows = num_rows;

  SLEQP_CALL(sleqp_calloc(&matrix->data, nnz_max));
  SLEQP_CALL(sleqp_calloc(&matrix->cols, num_cols + 1));
  SLEQP_CALL(sleqp_calloc(&matrix->rows, nnz_max));

  for(int i = 0; i < num_cols + 1; ++i)
  {
    matrix->cols[i] = 0;
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_matrix_reserve(SleqpSparseMatrix* matrix,
                                          int nnz)
{
  if(matrix->nnz_max >= nnz)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_realloc(&matrix->data, nnz));
  SLEQP_CALL(sleqp_realloc(&matrix->rows, nnz));

  matrix->nnz_max = nnz;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_matrix_resize(SleqpSparseMatrix* matrix,
                                         int num_rows,
                                         int num_cols)
{
  if(matrix->num_cols < num_cols)
  {
    SLEQP_CALL(sleqp_realloc(&matrix->cols, num_cols + 1));

    if(matrix->num_cols == 0)
    {
      matrix->cols[0] = 0;
    }

    for(int index = matrix->num_cols + 1; index < num_cols + 1; ++index)
    {
      matrix->cols[index] = matrix->cols[index - 1];
    }
  }
  else if(matrix->num_cols > num_cols)
  {
    matrix->nnz = matrix->cols[num_cols];
  }

  matrix->num_cols = num_cols;
  matrix->num_rows = num_rows;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_matrix_push(SleqpSparseMatrix* matrix,
                                       int row,
                                       int col,
                                       double value)
{
  assert(matrix->nnz < matrix->nnz_max);
  assert(row < matrix->num_rows);
  assert(col < matrix->num_cols);

  matrix->data[matrix->nnz] = value;
  matrix->rows[matrix->nnz] = row;
  matrix->cols[col + 1]++;
  matrix->nnz++;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_matrix_add_column(SleqpSparseMatrix* matrix,
                                             int col)
{
  assert(col < matrix->num_cols);

  matrix->cols[col + 1] = matrix->cols[col];

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_matrix_remove_column(SleqpSparseMatrix* matrix,
                                                int col)
{
  assert(col < matrix->num_cols);

  int nnz = matrix->cols[col + 1] - matrix->cols[col];

  matrix->cols[col + 1] = matrix->cols[col];
  matrix->nnz -= nnz;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_matrix_vector_product(SleqpSparseMatrix* matrix,
                                                 SleqpSparseVec* vector,
                                                 double* result)
{
  assert(matrix->num_cols == vector->dim);

  for(int index = 0; index < matrix->num_rows; ++index)
  {
    result[index] = 0.;
  }

  int k_vec = 0;

  while(k_vec < vector->nnz)
  {
    int col = vector->indices[k_vec];
    double factor = vector->data[k_vec];

    for(int entry = matrix->cols[col]; entry < matrix->cols[col + 1]; ++entry)
    {
      result[matrix->rows[entry]] += factor * matrix->data[entry];
    }

    ++k_vec;
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_matrix_trans_vector_product(SleqpSparseMatrix* matrix,
                                                       SleqpSparseVec* vector,
                                                       double eps,
                                                       SleqpSparseVec* result)
{
  assert(matrix->num_rows == vector->dim);
  assert(matrix->num_cols == result->dim);

  int col_size = 0;

  for(int col = 0; col < matrix->num_cols; ++col)
  {
    col_size += (matrix->cols[col + 1] - matrix->cols[col]) > 0;
  }

  SLEQP_CALL(sleqp_sparse_vector_clear(result));
  SLEQP_CALL(sleqp_sparse_vector_reserve(result, col_size));

  for(int col = 0; col < matrix->num_cols; ++col)
  {
    int k_vec = 0, k_mat = matrix->cols[col];

    double sum = 0.;

    while(k_vec < vector->nnz && k_mat < matrix->cols[col + 1])
    {
      int vec_idx = vector->indices[k_vec];
      int row_idx = matrix->rows[k_mat];

      if(vec_idx < row_idx)
      {
        ++k_vec;
      }
      else if(vec_idx > row_idx)
      {
        ++k_mat;
      }
      else
      {
        sum += vector->data[k_vec++] * matrix->data[k_mat++];
      }
    }

    if(!sleqp_zero(sum, eps))
    {
      SLEQP_CALL(sleqp_sparse_vector_push(result, col, sum));
    }

  }

  return SLEQP_OKAY;
}

double* sleqp_sparse_matrix_at(SleqpSparseMatrix* matrix,
                               int row,
                               int col)
{
  assert(row < matrix->num_rows);
  assert(col < matrix->num_cols);

  for(int index = matrix->cols[col]; index < matrix->cols[col + 1]; ++index)
  {
    if(matrix->rows[index] == row)
    {
      return matrix->data + index;
    }
  }

  return NULL;
}

SLEQP_RETCODE sleqp_sparse_matrix_clear(SleqpSparseMatrix* matrix)
{
  matrix->nnz = 0;

  for(int col = 0; col <= matrix->num_cols; ++col)
  {
    matrix->cols[col] = 0;
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_matrix_fprintf(SleqpSparseMatrix* matrix,
                                          FILE* output)
{
  fprintf(output,
          "Sparse matrix, dimension: %d x %d, entries: %d\n",
          matrix->num_rows,
          matrix->num_cols,
          matrix->nnz);

  int col = 0;

  for(int index = 0; index < matrix->nnz; ++index)
  {
    while(index >= matrix->cols[col + 1])
    {
      ++col;
    }

    assert(matrix->cols[col] <= index);
    assert(index < matrix->cols[col + 1]);

    fprintf(output, "(%d, %d) = %e\n",
            matrix->rows[index],
            col,
            matrix->data[index]);
  }

  return SLEQP_OKAY;
}

bool sleqp_sparse_matrix_valid(SleqpSparseMatrix* matrix)
{
  if(matrix->nnz > matrix->nnz_max)
  {
    return false;
  }

  if(matrix->num_cols < 0 || matrix->num_rows < 0)
  {
    return false;
  }

  if(matrix->nnz == 0)
  {
    return true;
  }

  for(int col = 0; col < matrix->num_cols; ++col)
  {
    if(matrix->cols[col] > matrix->cols[col + 1])
    {
      return false;
    }

    for(int index = matrix->cols[col]; index < matrix->cols[col + 1] - 1; ++index)
    {
      if(matrix->rows[index] >= matrix->rows[index + 1])
      {
        return false;
      }
    }

    for(int index = matrix->cols[col]; index < matrix->cols[col + 1]; ++index)
    {
      if(matrix->rows[index] < 0 ||
         matrix->rows[index] >= matrix->num_rows)
      {
        return false;
      }
    }
  }

  if(matrix->cols[matrix->num_cols] != matrix->nnz)
  {
    return false;
  }

  return true;
}

SLEQP_RETCODE sleqp_sparse_matrix_free(SleqpSparseMatrix** mstar)
{
  SleqpSparseMatrix* matrix = *mstar;

  sleqp_free(&matrix->rows);
  sleqp_free(&matrix->cols);
  sleqp_free(&matrix->data);

  sleqp_free(&matrix);

  *mstar = NULL;

  return SLEQP_OKAY;
}
