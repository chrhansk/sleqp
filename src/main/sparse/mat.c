#include "mat.h"

#include <math.h>

#include "cmp.h"
#include "error.h"
#include "fail.h"
#include "log.h"
#include "mem.h"

typedef struct SleqpMat
{
  int refcount;

  int num_rows;
  int num_cols;

  int nnz;
  int nnz_max;

  double* data;
  int* cols;
  int* rows;

} SleqpMat;

SLEQP_RETCODE
sleqp_mat_create(SleqpMat** mstar, int num_rows, int num_cols, int nnz_max)
{
  SLEQP_CALL(sleqp_malloc(mstar));

  SleqpMat* matrix = *mstar;

  *matrix = (SleqpMat){0};

  matrix->refcount = 1;

  matrix->nnz     = 0;
  matrix->nnz_max = nnz_max;

  matrix->num_cols = num_cols;
  matrix->num_rows = num_rows;

  SLEQP_CALL(sleqp_alloc_array(&matrix->data, nnz_max));
  SLEQP_CALL(sleqp_alloc_array(&matrix->cols, num_cols + 1));
  SLEQP_CALL(sleqp_alloc_array(&matrix->rows, nnz_max));

  for (int i = 0; i < num_cols + 1; ++i)
  {
    matrix->cols[i] = 0;
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_mat_reserve(SleqpMat* matrix, int nnz)
{
  assert(nnz >= 0);

  if (matrix->nnz_max >= nnz)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_realloc(&matrix->data, nnz));
  SLEQP_CALL(sleqp_realloc(&matrix->rows, nnz));

  matrix->nnz_max = nnz;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_mat_resize(SleqpMat* matrix, int num_rows, int num_cols)
{
  assert(num_rows >= 0);
  assert(num_cols >= 0);

  if (matrix->num_cols < num_cols)
  {
    SLEQP_CALL(sleqp_realloc(&matrix->cols, num_cols + 1));

    if (matrix->num_cols == 0)
    {
      matrix->cols[0] = 0;
    }

    for (int index = matrix->num_cols + 1; index < num_cols + 1; ++index)
    {
      matrix->cols[index] = matrix->cols[index - 1];
    }
  }
  else if (matrix->num_cols > num_cols)
  {
    matrix->nnz = matrix->cols[num_cols];
  }

  matrix->num_cols = num_cols;
  matrix->num_rows = num_rows;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_mat_scale(SleqpMat* matrix, double scale)
{
  if (scale == 0.)
  {
    SLEQP_CALL(sleqp_mat_clear(matrix));

    return SLEQP_OKAY;
  }

  for (int index = 0; index < matrix->nnz; ++index)
  {
    matrix->data[index] *= scale;
  }

  return SLEQP_OKAY;
}

int
sleqp_mat_num_cols(const SleqpMat* matrix)
{
  return matrix->num_cols;
}

int
sleqp_mat_num_rows(const SleqpMat* matrix)
{
  return matrix->num_rows;
}

int
sleqp_mat_nnz(const SleqpMat* matrix)
{
  return matrix->nnz;
}

int
sleqp_mat_nnz_max(const SleqpMat* matrix)
{
  return matrix->nnz_max;
}

SLEQP_RETCODE
sleqp_mat_set_nnz(SleqpMat* matrix, int nnz)
{
  matrix->nnz = nnz;
  return SLEQP_OKAY;
}

bool
sleqp_mat_is_quadratic(const SleqpMat* matrix)
{
  return matrix->num_rows == matrix->num_cols;
}

double*
sleqp_mat_data(const SleqpMat* matrix)
{
  return matrix->data;
}

int*
sleqp_mat_cols(const SleqpMat* matrix)
{
  return matrix->cols;
}

int*
sleqp_mat_rows(const SleqpMat* matrix)
{
  return matrix->rows;
}

SLEQP_RETCODE
sleqp_mat_push(SleqpMat* matrix, int row, int col, double value)
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

SLEQP_RETCODE
sleqp_mat_push_vec(SleqpMat* matrix, int col, SleqpVec* vec)
{
  assert(matrix->cols[col] == matrix->cols[col + 1]);
  assert(vec->dim == matrix->num_rows);

  assert((matrix->nnz_max - matrix->nnz) >= vec->nnz);

  for (int i = 0; i < vec->nnz; ++i)
  {
    matrix->data[matrix->nnz + i] = vec->data[i];
  }

  for (int i = 0; i < vec->nnz; ++i)
  {
    matrix->rows[matrix->nnz + i] = vec->indices[i];
  }

  matrix->cols[col + 1] += vec->nnz;
  matrix->nnz += vec->nnz;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_mat_push_col(SleqpMat* matrix, int col)
{
  assert(col < matrix->num_cols);

  matrix->cols[col + 1] = matrix->cols[col];

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_mat_pop_col(SleqpMat* matrix, int col)
{
  assert(col < matrix->num_cols);

  int nnz = matrix->cols[col + 1] - matrix->cols[col];

  matrix->cols[col + 1] = matrix->cols[col];
  matrix->nnz -= nnz;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_mat_vstack(const SleqpMat* first,
                 const SleqpMat* second,
                 SleqpMat* result)
{
  result->nnz = 0;

  assert(first->num_cols == second->num_cols);

  const int result_num_rows = first->num_rows + second->num_rows;
  const int result_num_cols = first->num_cols;

  SLEQP_CALL(sleqp_mat_resize(result, result_num_rows, result_num_cols));

  SLEQP_CALL(sleqp_mat_reserve(result, first->nnz + second->nnz));

  for (int col = 0; col < first->num_cols; ++col)
  {
    SLEQP_CALL(sleqp_mat_push_col(result, col));

    for (int k_first = first->cols[col]; k_first < first->cols[col + 1];
         ++k_first)
    {
      SLEQP_CALL(sleqp_mat_push(result,
                                first->rows[k_first],
                                col,
                                first->data[k_first]));
    }

    for (int k_second = second->cols[col]; k_second < second->cols[col + 1];
         ++k_second)
    {
      SLEQP_CALL(sleqp_mat_push(result,
                                first->num_rows + second->rows[k_second],
                                col,
                                second->data[k_second]));
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_mat_mult_vec(const SleqpMat* matrix,
                   const SleqpVec* vector,
                   double* result)
{
  assert(matrix->num_cols == vector->dim);

  for (int index = 0; index < matrix->num_rows; ++index)
  {
    result[index] = 0.;
  }

  int k_vec = 0;

  while (k_vec < vector->nnz)
  {
    int col       = vector->indices[k_vec];
    double factor = vector->data[k_vec];

    for (int entry = matrix->cols[col]; entry < matrix->cols[col + 1]; ++entry)
    {
      result[matrix->rows[entry]] += factor * matrix->data[entry];
    }

    ++k_vec;
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_mat_mult_vec_trans(const SleqpMat* matrix,
                         const SleqpVec* vector,
                         double eps,
                         SleqpVec* result)
{
  assert(matrix->num_rows == vector->dim);
  assert(matrix->num_cols == result->dim);

  int col_size = 0;

  for (int col = 0; col < matrix->num_cols; ++col)
  {
    col_size += (matrix->cols[col + 1] - matrix->cols[col]) > 0;
  }

  SLEQP_CALL(sleqp_vec_clear(result));
  SLEQP_CALL(sleqp_vec_reserve(result, col_size));

  for (int col = 0; col < matrix->num_cols; ++col)
  {
    int k_vec = 0, k_mat = matrix->cols[col];

    double sum = 0.;

    while (k_vec < vector->nnz && k_mat < matrix->cols[col + 1])
    {
      int vec_idx = vector->indices[k_vec];
      int row_idx = matrix->rows[k_mat];

      if (vec_idx < row_idx)
      {
        ++k_vec;
      }
      else if (vec_idx > row_idx)
      {
        ++k_mat;
      }
      else
      {
        sum += vector->data[k_vec++] * matrix->data[k_mat++];
      }
    }

    if (!sleqp_is_zero(sum, eps))
    {
      SLEQP_CALL(sleqp_vec_push(result, col, sum));
    }
  }

  return SLEQP_OKAY;
}

double*
sleqp_mat_at(SleqpMat* matrix, int row, int col)
{
  assert(row < matrix->num_rows);
  assert(col < matrix->num_cols);

  for (int index = matrix->cols[col]; index < matrix->cols[col + 1]; ++index)
  {
    if (matrix->rows[index] == row)
    {
      return matrix->data + index;
    }
  }

  return NULL;
}

double
sleqp_mat_value_at(SleqpMat* matrix, int row, int col)
{
  double* ptr = sleqp_mat_at(matrix, row, col);

  return ptr ? *ptr : 0.;
}

SLEQP_RETCODE
sleqp_mat_col(const SleqpMat* matrix, int col, SleqpVec* vec)
{
  assert(matrix->num_rows == vec->dim);
  assert(col >= 0);
  assert(col < matrix->num_cols);

  const int nnz = matrix->cols[col + 1] - matrix->cols[col];

  SLEQP_CALL(sleqp_vec_reserve(vec, nnz));
  SLEQP_CALL(sleqp_vec_clear(vec));

  for (int i = matrix->cols[col]; i < matrix->cols[col + 1]; ++i)
  {
    SLEQP_CALL(sleqp_vec_push(vec, matrix->rows[i], matrix->data[i]));
  }

  return SLEQP_OKAY;
}

bool
sleqp_mat_is_lower(const SleqpMat* matrix)
{
  for (int col = 0; col < matrix->num_cols; ++col)
  {
    for (int index = matrix->cols[col]; index < matrix->cols[col + 1]; ++index)
    {
      const int row = matrix->rows[index];

      if (row < col)
      {
        return false;
      }
    }
  }

  return true;
}

bool
sleqp_mat_eq(const SleqpMat* first, const SleqpMat* second, double eps)
{
  assert(first->num_rows == second->num_rows);
  assert(first->num_cols == second->num_cols);

  for (int col = 0; col < first->num_cols; ++col)
  {
    int first_index  = first->cols[col];
    int second_index = second->cols[col];

    bool first_valid  = true;
    bool second_valid = true;

    do
    {
      first_valid  = first_index < first->cols[col + 1];
      second_valid = second_index < second->cols[col + 1];

      bool both_valid = first_valid && second_valid;

      int first_row  = first_valid ? first->rows[first_index] : -1;
      int second_row = second_valid ? second->rows[second_index] : -1;

      double first_entry  = first_valid ? first->data[first_index] : -1;
      double second_entry = second_valid ? second->data[second_index] : -1;

      if (both_valid && (first_row == second_row))
      {
        if (!sleqp_is_eq(first_entry, second_entry, eps))
        {
          return false;
        }

        ++first_index;
        ++second_index;
      }
      else if (first_valid || (both_valid && (first_row < second_row)))
      {
        if (!sleqp_is_zero(first_entry, eps))
        {
          return false;
        }

        ++first_index;
      }
      else if (second_valid || (both_valid && (second_row < first_row)))
      {
        if (!sleqp_is_zero(second_entry, eps))
        {
          return false;
        }

        ++second_index;
      }
    } while (first_valid || second_valid);
  }

  return true;
}

SLEQP_RETCODE
sleqp_mat_remove_rows(const SleqpMat* source,
                      SleqpMat* target,
                      const int* row_indices,
                      int num_row_entries)
{
  SLEQP_CALL(sleqp_mat_reserve(target, sleqp_mat_nnz(source)));

  SLEQP_CALL(sleqp_mat_clear(target));

  assert(sleqp_mat_num_cols(source) == sleqp_mat_num_cols(target));

  assert(sleqp_mat_num_rows(source)
         == sleqp_mat_num_rows(target) + num_row_entries);

  assert(row_indices || (num_row_entries == 0));

  const int num_cols = sleqp_mat_num_cols(source);

  double* source_data = sleqp_mat_data(source);
  int* source_rows    = sleqp_mat_rows(source);
  int* source_cols    = sleqp_mat_cols(source);

  for (int col = 0; col < num_cols; ++col)
  {
    SLEQP_CALL(sleqp_mat_push_col(target, col));

    int row_offset = 0;

    for (int k = source_cols[col]; k < source_cols[col + 1]; ++k)
    {
      const int row = source_rows[k];

      while (row_offset < num_row_entries && row_indices[row_offset] < row)
      {
        ++row_offset;
      }

      if (row_offset < num_row_entries && row_indices[row_offset] == row)
      {
        continue;
      }

      SLEQP_CALL(sleqp_mat_push(target, row - row_offset, col, source_data[k]));
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_mat_remove_entries(const SleqpMat* source,
                         SleqpMat* target,
                         const int* col_indices,
                         int num_col_entries,
                         const int* row_indices,
                         int num_row_entries)
{
  SLEQP_CALL(sleqp_mat_reserve(target, sleqp_mat_nnz(source)));

  SLEQP_CALL(sleqp_mat_clear(target));

  assert(sleqp_mat_num_cols(source)
         == sleqp_mat_num_cols(target) + num_col_entries);

  assert(sleqp_mat_num_rows(source)
         == sleqp_mat_num_rows(target) + num_row_entries);

  assert(row_indices || (num_row_entries == 0));
  assert(col_indices || (num_col_entries == 0));

  const int num_cols = sleqp_mat_num_cols(source);

  double* source_data = sleqp_mat_data(source);
  int* source_rows    = sleqp_mat_rows(source);
  int* source_cols    = sleqp_mat_cols(source);

  int col_offset = 0;

  for (int col = 0; col < num_cols; ++col)
  {
    if (col_offset < num_col_entries && col_indices[col_offset] <= col)
    {
      ++col_offset;
      continue;
    }

    SLEQP_CALL(sleqp_mat_push_col(target, col - col_offset));

    int row_offset = 0;

    for (int k = source_cols[col]; k < source_cols[col + 1]; ++k)
    {
      const int row = source_rows[k];

      while (row_offset < num_row_entries && row_indices[row_offset] < row)
      {
        ++row_offset;
      }

      if (row_offset < num_row_entries && row_indices[row_offset] == row)
      {
        continue;
      }

      SLEQP_CALL(sleqp_mat_push(target,
                                row - row_offset,
                                col - col_offset,
                                source_data[k]));
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_mat_remove_cols(const SleqpMat* source,
                      SleqpMat* target,
                      const int* col_indices,
                      int num_col_entries)
{
  SLEQP_CALL(sleqp_mat_reserve(target, sleqp_mat_nnz(source)));

  SLEQP_CALL(sleqp_mat_clear(target));

  assert(sleqp_mat_num_cols(source)
         == sleqp_mat_num_cols(target) + num_col_entries);

  assert(sleqp_mat_num_rows(source) == sleqp_mat_num_rows(target));

  assert(col_indices || (num_col_entries == 0));

  const int num_cols = sleqp_mat_num_cols(source);

  double* source_data = sleqp_mat_data(source);
  int* source_rows    = sleqp_mat_rows(source);
  int* source_cols    = sleqp_mat_cols(source);

  int col_offset = 0;

  for (int col = 0; col < num_cols; ++col)
  {
    if (col_offset < num_col_entries && col_indices[col_offset] <= col)
    {
      ++col_offset;
      continue;
    }

    SLEQP_CALL(sleqp_mat_push_col(target, col - col_offset));

    for (int k = source_cols[col]; k < source_cols[col + 1]; ++k)
    {
      const int row = source_rows[k];

      SLEQP_CALL(sleqp_mat_push(target, row, col - col_offset, source_data[k]));
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_mat_clear(SleqpMat* matrix)
{
  matrix->nnz = 0;

  for (int col = 0; col <= matrix->num_cols; ++col)
  {
    matrix->cols[col] = 0;
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_mat_fprintf(const SleqpMat* matrix, FILE* output)
{
  fprintf(output,
          "Sparse matrix, dimension: %d x %d, entries: %d\n",
          matrix->num_rows,
          matrix->num_cols,
          matrix->nnz);

  int col = 0;

  for (int index = 0; index < matrix->nnz; ++index)
  {
    while (index >= matrix->cols[col + 1])
    {
      ++col;
    }

    assert(matrix->cols[col] <= index);
    assert(index < matrix->cols[col + 1]);

    fprintf(output,
            "(%d, %d) = %e\n",
            matrix->rows[index],
            col,
            matrix->data[index]);
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_mat_dump(const SleqpMat* matrix, FILE* output)
{
  fprintf(output, "%%%%MatrixMarket matrix coordinate real general\n");
  fprintf(output,
          "%d %d %d\n",
          matrix->num_rows,
          matrix->num_cols,
          matrix->nnz);

  int col = 0;

  for (int index = 0; index < matrix->nnz; ++index)
  {
    while (index >= matrix->cols[col + 1])
    {
      ++col;
    }

    fprintf(output,
            "%d %d %f\n",
            matrix->rows[index] + 1,
            col + 1,
            matrix->data[index]);
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_mat_dump_to_file(const SleqpMat* matrix, const char* name)
{
  FILE* output = fopen(name, "w");

  if (!output)
  {
    sleqp_raise(SLEQP_ILLEGAL_ARGUMENT,
                "Failed to open output file '%s'",
                name);
  }

  SLEQP_CALL(sleqp_mat_dump(matrix, output));

  fclose(output);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_mat_copy(const SleqpMat* source, SleqpMat* target)
{
  assert(source->num_cols == target->num_cols);
  assert(source->num_rows == target->num_rows);
  assert(sleqp_mat_is_valid(source));

  SLEQP_CALL(sleqp_mat_clear(target));

  SLEQP_CALL(sleqp_mat_reserve(target, source->nnz));

  for (int i = 0; i < source->nnz; ++i)
  {
    target->data[i] = source->data[i];
    target->rows[i] = source->rows[i];
  }

  for (int i = 0; i <= target->num_cols; ++i)
  {
    target->cols[i] = source->cols[i];
  }

  target->nnz = source->nnz;

  assert(sleqp_mat_is_valid(target));

  return SLEQP_OKAY;
}

bool
sleqp_mat_is_valid(const SleqpMat* matrix)
{
  if (matrix->nnz > matrix->nnz_max)
  {
    return false;
  }

  if (matrix->num_cols < 0 || matrix->num_rows < 0)
  {
    return false;
  }

  if (matrix->nnz == 0)
  {
    return true;
  }

  for (int col = 0; col < matrix->num_cols; ++col)
  {
    if (matrix->cols[col] > matrix->cols[col + 1])
    {
      return false;
    }

    for (int index = matrix->cols[col]; index < matrix->cols[col + 1] - 1;
         ++index)
    {
      if (matrix->rows[index] >= matrix->rows[index + 1])
      {
        return false;
      }
    }

    for (int index = matrix->cols[col]; index < matrix->cols[col + 1]; ++index)
    {
      if (matrix->rows[index] < 0 || matrix->rows[index] >= matrix->num_rows)
      {
        return false;
      }
    }
  }

  if (matrix->cols[matrix->num_cols] != matrix->nnz)
  {
    return false;
  }

  return true;
}

static SLEQP_RETCODE
compute_row_sums(const SleqpMat* matrix, int* row_sums)
{
  const int matrix_num_rows = sleqp_mat_num_rows(matrix);
  const int matrix_nnz      = sleqp_mat_nnz(matrix);

  const int* matrix_rows = sleqp_mat_rows(matrix);

  for (int i = 0; i < matrix_num_rows; ++i)
  {
    row_sums[i] = 0;
  }

  for (int k = 0; k < matrix_nnz; ++k)
  {
    ++(row_sums[matrix_rows[k]]);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
fill_transpose(const SleqpMat* source, SleqpMat* target, int* row_sums)
{
  const int nnz = sleqp_mat_nnz(source);

  SLEQP_CALL(sleqp_mat_clear(target));
  SLEQP_CALL(sleqp_mat_reserve(target, nnz));

  {
    const int num_rows = sleqp_mat_num_rows(source);
    const int num_cols = sleqp_mat_num_cols(source);
    SLEQP_CALL(sleqp_mat_resize(target, num_cols, num_rows));
  }

  {
    int* matrix_cols   = sleqp_mat_cols(target);
    const int num_cols = sleqp_mat_num_cols(target);

    matrix_cols[0] = 0;
    int offset     = 0;

    for (int j = 0; j < num_cols; ++j)
    {
      offset += row_sums[j];
      matrix_cols[j + 1] = offset;
      row_sums[j]        = 0;
    }
  }

  // Store elements
  {
    const int* source_rows    = sleqp_mat_rows(source);
    const int* source_cols    = sleqp_mat_cols(source);
    const double* source_data = sleqp_mat_data(source);
    const int source_num_cols = sleqp_mat_num_cols(source);

    int* target_rows    = sleqp_mat_rows(target);
    int* target_cols    = sleqp_mat_cols(target);
    double* target_data = sleqp_mat_data(target);

    for (int col = 0; col < source_num_cols; ++col)
    {
      for (int k = source_cols[col]; k < source_cols[col + 1]; ++k)
      {
        const double value = source_data[k];
        const double row   = source_rows[k];

        const int target_col = row;
        const int target_row = col;
        const int target_pos = target_cols[target_col] + row_sums[target_col];
        ++(row_sums[target_col]);

        target_data[target_pos] = value;
        target_rows[target_pos] = target_row;
      }
    }
  }

  // Update nnz
  SLEQP_CALL(sleqp_mat_set_nnz(target, nnz));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_mat_trans(const SleqpMat* source, SleqpMat* target, int* row_cache)
{
  SLEQP_CALL(compute_row_sums(source, row_cache));
  SLEQP_CALL(fill_transpose(source, target, row_cache));

  assert(sleqp_mat_is_valid(target));

  return SLEQP_OKAY;
}

bool
sleqp_mat_is_finite(const SleqpMat* matrix)
{
  for (int index = 0; index < matrix->nnz; ++index)
  {
    if (!sleqp_is_finite(matrix->data[index]))
    {
      return false;
    }
  }

  return true;
}

static SLEQP_RETCODE
sparse_matrix_free(SleqpMat** mstar)
{
  SleqpMat* matrix = *mstar;

  if (!matrix)
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

SLEQP_RETCODE
sleqp_mat_capture(SleqpMat* matrix)
{
  ++matrix->refcount;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_mat_release(SleqpMat** star)
{
  SleqpMat* matrix = *star;

  if (!matrix)
  {
    return SLEQP_OKAY;
  }

  if (--matrix->refcount == 0)
  {
    SLEQP_CALL(sparse_matrix_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
