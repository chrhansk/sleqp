#include "sleqp_sparse.h"

#include <assert.h>

#include "sleqp_cmp.h"
#include "sleqp_mem.h"

SLEQP_RETCODE sleqp_sparse_vector_create(SleqpSparseVec** vstar,
                                         size_t dim,
                                         size_t nnz_max)
{
  assert(nnz_max <= dim);

  SLEQP_CALL(sleqp_malloc(vstar));

  SleqpSparseVec *vec = *vstar;

  vec->nnz = 0;
  vec->dim = dim;
  vec->nnz_max = nnz_max;

  SLEQP_CALL(sleqp_calloc(&vec->data, nnz_max));
  SLEQP_CALL(sleqp_calloc(&vec->indices, nnz_max));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_vector_push(SleqpSparseVec* vec,
                                       size_t idx,
                                       double value)
{
  assert(idx < vec->nnz_max);

  if(vec->nnz > 0)
  {
    assert(idx > vec->indices[vec->nnz - 1]);
  }

  vec->data[vec->nnz] = value;
  vec->indices[vec->nnz] = idx;

  ++(vec->nnz);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_vector_from_raw(SleqpSparseVec* vec,
                                           double* values,
                                           size_t dim)
{
  size_t nnz = 0;

  for(size_t i = 0; i < dim;++i)
  {
    if(!sleqp_zero(values[i]))
    {
      ++nnz;
    }
  }

  vec->dim = dim;
  sleqp_sparse_vector_reserve(vec, nnz);

  for(size_t i = 0; i < dim;++i)
  {
    double v = values[i];

    if(!sleqp_zero(v))
    {
      sleqp_sparse_vector_push(vec, i, v);
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_vector_reserve(SleqpSparseVec* vec,
                                          size_t nnz_max)
{
  if(vec->nnz_max >= nnz_max)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_realloc(&vec->data, nnz_max));
  SLEQP_CALL(sleqp_realloc(&vec->indices, nnz_max));

  vec->nnz_max = nnz_max;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_vector_dot(SleqpSparseVec* first,
                                      SleqpSparseVec* second,
                                      double* product)
{
  assert(first->dim == second->dim);

  *product = 0.;

  size_t k_first = 0, k_second = 0;

  while(k_first < first->nnz && k_second < second->nnz)
  {
    size_t i_first = first->indices[k_first];
    size_t i_second = second->indices[k_second];

    if(i_first == i_second)
    {
      *product += first->data[k_first] * second->data[k_second];
      ++k_first;
      ++k_second;
    }
    else if(i_first < i_second)
    {
      ++k_first;
    }
    else
    {
      ++k_second;
    }

  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_vector_scale(SleqpSparseVec* vector,
                                        double factor)
{
  for(size_t k = 0; k < vector->nnz; ++k)
  {
    vector->data[k] *= factor;
  }
}

SLEQP_RETCODE sleqp_sparse_vector_dense_dot(SleqpSparseVec* first,
                                            double* second,
                                            double* product)
{
  *product = 0.;

  for(size_t k_first = 0;k_first < first->nnz; ++k_first)
  {
    size_t i_first = first->indices[k_first];

    *product += first->data[k_first] * second[i_first];

    ++k_first;
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_vector_add(SleqpSparseVec* first,
                                      SleqpSparseVec* second,
                                      double first_factor,
                                      double second_factor,
                                      SleqpSparseVec* result)
{
  assert(first->dim == second->dim);

  result->dim = first->dim;
  result->nnz = 0;

  size_t k_first = 0, k_second = 0;

  while(k_first < first->nnz || k_second < second->nnz)
  {
    SLEQP_Bool valid_first = (k_first < first->nnz);
    SLEQP_Bool valid_second = (k_second < second->nnz);

    size_t i_first = valid_first ? first->indices[k_first] : first->dim + 1;
    size_t i_second = valid_second ? second->indices[k_second] : first->dim + 1;

    size_t i_combined = SLEQP_MIN(i_first, i_second);

    double value = 0.;

    if(i_first == i_combined)
    {
      value += first_factor * first->data[k_first];
      ++k_first;
    }

    if(i_second == i_combined)
    {
      value += second_factor * second->data[k_second];
      ++k_second;
    }

    SLEQP_CALL(sleqp_sparse_vector_push(result,
                                        i_combined,
                                        value));
  }

  return SLEQP_OKAY;
}

double* sleqp_sparse_vector_at(SleqpSparseVec* vec,
                               size_t index)
{
  assert(index < vec->dim);

  for(int i = 0; i < vec->nnz; ++i)
  {
    if(vec->indices[i] == index)
    {
      return vec->data + i;
    }
  }

  return NULL;
}

SLEQP_RETCODE sleqp_sparse_vector_clip(SleqpSparseVec* x,
                                       SleqpSparseVec* lb,
                                       SleqpSparseVec* ub,
                                       SleqpSparseVec** xstar)
{
  const size_t dim = x->dim;

  assert(lb->dim == dim);
  assert(ub->dim == dim);

  sleqp_sparse_vector_create(xstar,
                             dim,
                             SLEQP_MIN(x->nnz + lb->nnz, dim));

  size_t k_x = 0, k_lb = 0, k_ub = 0;

  SleqpSparseVec* xclip = *xstar;

  while(k_x < x->nnz || k_lb < lb->nnz || k_ub < ub->nnz)
  {
    double x_val = 0;

    SLEQP_Bool valid_x = (k_x < x->nnz);
    SLEQP_Bool valid_lb = (k_lb < lb->nnz);
    SLEQP_Bool valid_ub = (k_ub < ub->nnz);

    size_t idx = valid_x ? x->indices[k_x] : dim + 1;
    idx = SLEQP_MIN(idx, valid_lb ? lb->indices[k_lb] : dim + 1);
    idx = SLEQP_MIN(idx, valid_ub ? ub->indices[k_ub] : dim + 1);

    if(valid_x && idx == x->indices[k_x])
    {
      x_val = x->data[k_x];
    }
    else
    {
      x_val = 0.;
    }

    if(valid_lb && idx == lb->indices[k_lb])
    {
      x_val = SLEQP_MAX(x_val, lb->data[k_lb]);
    }

    if(valid_ub && idx == ub->indices[k_ub])
    {
      x_val = SLEQP_MIN(x_val, ub->data[k_ub]);
    }

    if(!sleqp_zero(x_val))
    {
      SLEQP_CALL(sleqp_sparse_vector_push(xclip,
                                          idx,
                                          x_val));
    }

    if(valid_lb && idx == lb->indices[k_lb])
    {
      ++k_lb;
    }

    if(valid_ub && idx == ub->indices[k_ub])
    {
      ++k_ub;
    }

    if(valid_x && idx == x->indices[k_x])
    {
      ++k_x;
    }

  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_vector_fprintf(SleqpSparseVec* vec,
                                          FILE* output)
{
  fprintf(output,
          "Sparse vector, dimension: %ld, entries: %ld\n",
          vec->dim,
          vec->nnz);

  for(size_t index = 0; index < vec->nnz; ++index)
  {
    fprintf(output, "(%ld) = %f\n",
            vec->indices[index],
            vec->data[index]);
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_vector_free(SleqpSparseVec** vstar)
{
  SleqpSparseVec *vec = *vstar;

  sleqp_free(&vec->indices);
  sleqp_free(&vec->data);

  sleqp_free(&vec);

  *vstar = NULL;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_matrix_create(SleqpSparseMatrix** mstar,
                                         size_t num_rows,
                                         size_t num_cols,
                                         size_t nnz_max)
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
                                          size_t nnz)
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
                                         size_t num_rows,
                                         size_t num_cols)
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
    for(int column = matrix->num_cols - 1; column >= num_cols; --column)
    {
      SLEQP_CALL(sleqp_sparse_matrix_remove_column(matrix, column));
    }
  }

  matrix->num_cols = num_cols;
  matrix->num_rows = num_rows;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_matrix_push(SleqpSparseMatrix* matrix,
                                       size_t row,
                                       size_t col,
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
                                             size_t col)
{
  assert(col < matrix->num_cols);

  matrix->cols[col + 1] = matrix->cols[col];

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_matrix_remove_column(SleqpSparseMatrix* matrix,
                                                size_t col)
{
  assert(col < matrix->num_cols);

  size_t nnz = matrix->cols[col + 1] - matrix->cols[col];

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

  size_t k_vec = 0.;

  while(k_vec < vector->nnz)
  {
    int col = vector->indices[k_vec];
    double factor = vector->data[k_vec];

    for(int entry = matrix->cols[col]; entry< matrix->cols[col + 1]; ++entry)
    {
      result[matrix->rows[entry]] += factor * matrix->data[entry];
    }
  }

  return SLEQP_OKAY;
}

double* sleqp_sparse_matrix_at(SleqpSparseMatrix* matrix,
                               size_t row,
                               size_t col)
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

    fprintf(output, "(%d, %d) = %f\n",
            matrix->rows[index],
            col,
            matrix->data[index]);
  }

  return SLEQP_OKAY;
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
