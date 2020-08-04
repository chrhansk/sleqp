#include "sleqp_sparse_matrix.h"

#include <assert.h>

#include "sleqp_cmp.h"
#include "sleqp_mem.h"

typedef struct SleqpSparseMatrix
{
  int refcount;

  int num_rows;
  int num_cols;

  int nnz;
  int nnz_max;

  double* data;
  int* cols;
  int* rows;

} SleqpSparseMatrix;

SLEQP_RETCODE sleqp_sparse_matrix_create(SleqpSparseMatrix** mstar,
                                         int num_rows,
                                         int num_cols,
                                         int nnz_max)
{
  SLEQP_CALL(sleqp_malloc(mstar));

  SleqpSparseMatrix* matrix = *mstar;

  *matrix = (SleqpSparseMatrix){0};

  matrix->refcount = 1;

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

int sleqp_sparse_matrix_get_num_cols(const SleqpSparseMatrix* matrix)
{
  return matrix->num_cols;
}

int sleqp_sparse_matrix_get_num_rows(const SleqpSparseMatrix* matrix)
{
  return matrix->num_rows;
}

int sleqp_sparse_matrix_get_nnz(const SleqpSparseMatrix* matrix)
{
  return matrix->nnz;
}

int sleqp_sparse_matrix_get_nnz_max(const SleqpSparseMatrix* matrix)
{
  return matrix->nnz_max;
}

SLEQP_RETCODE sleqp_sparse_matrix_set_nnz(SleqpSparseMatrix* matrix,
                                          int nnz)
{
  matrix->nnz = nnz;
  return SLEQP_OKAY;
}

bool sleqp_sparse_matrix_is_quadratic(const SleqpSparseMatrix* matrix)
{
  return matrix->num_rows == matrix->num_cols;
}

double* sleqp_sparse_matrix_get_data(SleqpSparseMatrix* matrix)
{
  return matrix->data;
}

int* sleqp_sparse_matrix_get_cols(SleqpSparseMatrix* matrix)
{
  return matrix->cols;
}

int* sleqp_sparse_matrix_get_rows(SleqpSparseMatrix* matrix)
{
  return matrix->rows;
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

SLEQP_RETCODE sleqp_sparse_matrix_push_column(SleqpSparseMatrix* matrix,
                                              int col)
{
  assert(col < matrix->num_cols);

  matrix->cols[col + 1] = matrix->cols[col];

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_matrix_pop_column(SleqpSparseMatrix* matrix,
                                             int col)
{
  assert(col < matrix->num_cols);

  int nnz = matrix->cols[col + 1] - matrix->cols[col];

  matrix->cols[col + 1] = matrix->cols[col];
  matrix->nnz -= nnz;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_matrix_vector_product(const SleqpSparseMatrix* matrix,
                                                 const SleqpSparseVec* vector,
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

SLEQP_RETCODE sleqp_sparse_matrix_trans_vector_product(const SleqpSparseMatrix* matrix,
                                                       const SleqpSparseVec* vector,
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

bool sleqp_sparse_matrix_eq(const SleqpSparseMatrix* first,
                            const SleqpSparseMatrix* second,
                            double eps)
{
  assert(first->num_rows == second->num_rows);
  assert(first->num_cols == second->num_cols);

  for(int col = 0; col < first->num_cols; ++col)
  {
    int first_index = first->cols[col];
    int second_index = second->cols[col];

    bool first_valid = true;
    bool second_valid = true;

    do
    {
      first_valid = first_index < first->cols[col + 1];
      second_valid = second_index < second->cols[col + 1];

      bool both_valid = first_valid && second_valid;

      int first_row = first_valid ? first->rows[first_index] : -1;
      int second_row = second_valid ? second->rows[second_index] : -1;

      double first_entry = first_valid ? first->data[first_index] : -1;
      double second_entry = second_valid ? second->data[second_index] : -1;

      if(both_valid && (first_row == second_row))
      {
        if(!sleqp_eq(first_entry, second_entry, eps))
        {
          return false;
        }

        ++first_index;
        ++second_index;
      }
      else if(first_valid || (both_valid && (first_row < second_row)))
      {
        if(!sleqp_zero(first_entry, eps))
        {
          return false;
        }

        ++first_index;
      }
      else if(second_valid || (both_valid && (second_row < first_row)))
      {
        if(!sleqp_zero(second_entry, eps))
        {
          return false;
        }

        ++second_index;
      }
    }
    while(first_valid || second_valid);

  }

  return true;
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

SLEQP_RETCODE sleqp_sparse_matrix_fprintf(const SleqpSparseMatrix* matrix,
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

SLEQP_RETCODE sleqp_sparse_matrix_dump(const SleqpSparseMatrix* matrix,
                                       FILE* output)
{
  fprintf(output, "%%%%MatrixMarket matrix coordinate real general\n");
  fprintf(output, "%d %d %d\n", matrix->num_rows, matrix->num_cols, matrix->nnz);

  int col = 0;

  for(int index = 0; index < matrix->nnz; ++index)
  {
    while(index >= matrix->cols[col + 1])
    {
      ++col;
    }

    fprintf(output, "%d %d %f\n",
            matrix->rows[index] + 1,
            col + 1,
            matrix->data[index]);
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_matrix_copy(const SleqpSparseMatrix* source,
                                       SleqpSparseMatrix* target)
{
  assert(source->num_cols == target->num_cols);
  assert(source->num_rows == target->num_rows);
  assert(sleqp_sparse_matrix_valid(source));

  SLEQP_CALL(sleqp_sparse_matrix_clear(target));

  SLEQP_CALL(sleqp_sparse_matrix_reserve(target,
                                         source->nnz));

  for(int i = 0; i < source->nnz;++i)
  {
    target->data[i] = source->data[i];
    target->rows[i] = source->rows[i];
  }

  for(int i = 0; i <= target->num_cols; ++i)
  {
    target->cols[i] = source->cols[i];
  }

  target->nnz = source->nnz;

  assert(sleqp_sparse_matrix_valid(target));

  return SLEQP_OKAY;
}

bool sleqp_sparse_matrix_valid(const SleqpSparseMatrix* matrix)
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

static SLEQP_RETCODE sparse_matrix_free(SleqpSparseMatrix** mstar)
{
  SleqpSparseMatrix* matrix = *mstar;

  if(!matrix)
  {
    return SLEQP_OKAY;
  }

  sleqp_free(&matrix->rows);
  sleqp_free(&matrix->cols);
  sleqp_free(&matrix->data);

  sleqp_free(&matrix);

  *mstar = NULL;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_matrix_capture(SleqpSparseMatrix* matrix)
{
  ++matrix->refcount;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_matrix_release(SleqpSparseMatrix** star)
{
  SleqpSparseMatrix* matrix = *star;

  if(!matrix)
  {
    return SLEQP_OKAY;
  }

  if(--matrix->refcount == 0)
  {
    SLEQP_CALL(sparse_matrix_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
