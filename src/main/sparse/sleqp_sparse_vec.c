#include "sleqp_sparse_vec.h"

#include <assert.h>

#include "sleqp_cmp.h"
#include "sleqp_mem.h"

SLEQP_RETCODE sleqp_sparse_vector_create(SleqpSparseVec** vstar,
                                         int dim,
                                         int nnz_max)
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
                                       int idx,
                                       double value)
{
  assert(vec->nnz < vec->nnz_max);

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
                                           int dim)
{
  int nnz = 0;

  for(int i = 0; i < dim;++i)
  {
    if(!sleqp_zero(values[i]))
    {
      ++nnz;
    }
  }

  vec->nnz = 0;
  vec->dim = dim;

  SLEQP_CALL(sleqp_sparse_vector_reserve(vec, nnz));

  for(int i = 0; i < dim;++i)
  {
    double v = values[i];

    if(!sleqp_zero(v))
    {
      sleqp_sparse_vector_push(vec, i, v);
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_vector_to_raw(SleqpSparseVec* vec,
                                         double* values)
{
  for(int i = 0; i < vec->dim; ++i)
  {
    values[i] = 0.;
  }

  for(int k = 0; k < vec->nnz; ++k)
  {
    values[vec->indices[k]] = vec->data[k];
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_vector_copy(SleqpSparseVec* source,
                                       SleqpSparseVec* target)
{
  assert(source->dim == target->dim);

  SLEQP_CALL(sleqp_sparse_vector_reserve(target, source->nnz));

  target->nnz = 0;

  for(int k = 0; k < source->nnz; ++k)
  {
    SLEQP_CALL(sleqp_sparse_vector_push(target,
                                        source->indices[k],
                                        source->data[k]));
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_vector_clear(SleqpSparseVec* vec)
{
  vec->nnz = 0;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_vector_reserve(SleqpSparseVec* vec,
                                          int nnz_max)
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

SLEQP_RETCODE sleqp_sparse_vector_resize(SleqpSparseVec* vec,
                                         int dim)
{
  if(dim < vec->dim)
  {
    while(vec->nnz > 0 && vec->indices[vec->nnz - 1] >= dim)
    {
      --vec->nnz;
    }
  }

  vec->dim = dim;

  return SLEQP_OKAY;
}

bool sleqp_sparse_vector_eq(SleqpSparseVec* first,
                            SleqpSparseVec* second)
{
  assert(first->dim == second->dim);

  int k_first = 0, k_second = 0;

  while(k_first < first->nnz || k_second < second->nnz)
  {
    bool valid_first = (k_first < first->nnz);
    bool valid_second = (k_second < second->nnz);

    double first_value = valid_first ? first->data[k_first] : 0.;
    double second_value = valid_second ? second->data[k_second] : 0.;

    if(!sleqp_eq(first_value, second_value))
    {
      return false;
    }

    if(valid_first)
    {
      ++k_first;
    }

    if(valid_second)
    {
      ++k_second;
    }

  }

  return true;
}

SLEQP_RETCODE sleqp_sparse_vector_dot(SleqpSparseVec* first,
                                      SleqpSparseVec* second,
                                      double* product)
{
  assert(first->dim == second->dim);

  *product = 0.;

  int k_first = 0, k_second = 0;

  while(k_first < first->nnz && k_second < second->nnz)
  {
    int i_first = first->indices[k_first];
    int i_second = second->indices[k_second];

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
  for(int k = 0; k < vector->nnz; ++k)
  {
    vector->data[k] *= factor;
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_vector_dense_dot(SleqpSparseVec* first,
                                            double* second,
                                            double* product)
{
  *product = 0.;

  for(int k_first = 0;k_first < first->nnz; ++k_first)
  {
    int i_first = first->indices[k_first];

    *product += first->data[k_first] * second[i_first];

    ++k_first;
  }

  return SLEQP_OKAY;
}

double sleqp_sparse_vector_normsq(SleqpSparseVec* vec)
{
  double normsq = 0.;

  for(int k = 0; k < vec->nnz; ++k)
  {
    double value = vec->data[k];

    normsq += value*value;
  }

  return normsq;
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

  SLEQP_CALL(sleqp_sparse_vector_reserve(result, first->nnz + second->nnz));

  int k_first = 0, k_second = 0;

  while(k_first < first->nnz || k_second < second->nnz)
  {
    bool valid_first = (k_first < first->nnz);
    bool valid_second = (k_second < second->nnz);

    int i_first = valid_first ? first->indices[k_first] : first->dim + 1;
    int i_second = valid_second ? second->indices[k_second] : first->dim + 1;

    int i_combined = SLEQP_MIN(i_first, i_second);

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
                               int index)
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
  const int dim = x->dim;

  assert(lb->dim == dim);
  assert(ub->dim == dim);

  sleqp_sparse_vector_create(xstar,
                             dim,
                             SLEQP_MIN(x->nnz + lb->nnz, dim));

  int k_x = 0, k_lb = 0, k_ub = 0;

  SleqpSparseVec* xclip = *xstar;

  while(k_x < x->nnz || k_lb < lb->nnz || k_ub < ub->nnz)
  {
    double x_val = 0;

    bool valid_x = (k_x < x->nnz);
    bool valid_lb = (k_lb < lb->nnz);
    bool valid_ub = (k_ub < ub->nnz);

    int idx = valid_x ? x->indices[k_x] : dim + 1;
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
          "Sparse vector, dimension: %d, entries: %d\n",
          vec->dim,
          vec->nnz);

  for(int index = 0; index < vec->nnz; ++index)
  {
    fprintf(output, "(%d) = %f\n",
            vec->indices[index],
            vec->data[index]);
  }

  return SLEQP_OKAY;
}

bool sleqp_sparse_vector_valid(SleqpSparseVec* vec)
{
  if(vec->nnz > vec->nnz_max || vec->nnz < 0)
  {
    return false;
  }

  if(vec->nnz == 0)
  {
    return true;
  }

  for(int k = 0; k < vec->nnz; ++k)
  {
    if(vec->indices[k] < 0)
    {
      return false;
    }
  }

  for(int k = 0; k < vec->nnz - 1; ++k)
  {
    if(vec->indices[k] >= vec->indices[k + 1])
    {
      return false;
    }
  }

  if(vec->indices[vec->nnz - 1] >= vec->dim)
  {
    return false;
  }

  return true;
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
