#include "sparse_vec.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "cmp.h"
#include "error.h"
#include "log.h"
#include "mem.h"

SLEQP_RETCODE
sleqp_sparse_vector_create(SleqpSparseVec** vstar, int dim, int nnz_max)
{
  assert(dim >= 0);
  assert(nnz_max >= 0);

  SLEQP_CALL(sleqp_malloc(vstar));

  SleqpSparseVec* vec = *vstar;

  vec->nnz     = 0;
  vec->dim     = dim;
  vec->nnz_max = nnz_max;

  SLEQP_CALL(sleqp_alloc_array(&vec->data, nnz_max));
  SLEQP_CALL(sleqp_alloc_array(&vec->indices, nnz_max));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_sparse_vector_create_empty(SleqpSparseVec** vec, int dim)
{
  SLEQP_CALL(sleqp_sparse_vector_create(vec, dim, 0));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_sparse_vector_create_full(SleqpSparseVec** vec, int dim)
{
  SLEQP_CALL(sleqp_sparse_vector_create(vec, dim, dim));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_sparse_vector_push(SleqpSparseVec* vec, int idx, double value)
{
  assert(vec->nnz < vec->nnz_max);

  // assert(!(isinf(value) || isnan(value)));

  assert(idx < vec->dim);
  assert(idx >= 0);

  if (vec->nnz > 0)
  {
    assert(idx > vec->indices[vec->nnz - 1]);
  }

  vec->data[vec->nnz]    = value;
  vec->indices[vec->nnz] = idx;

  ++(vec->nnz);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_sparse_vector_from_raw(SleqpSparseVec* vec,
                             double* values,
                             int dim,
                             double zero_eps)
{
  int nnz = 0;

  for (int i = 0; i < dim; ++i)
  {
    if (!sleqp_is_zero(values[i], zero_eps))
    {
      ++nnz;
    }
  }

  SLEQP_CALL(sleqp_sparse_vector_clear(vec));
  SLEQP_CALL(sleqp_sparse_vector_resize(vec, dim));

  SLEQP_CALL(sleqp_sparse_vector_reserve(vec, nnz));

  for (int i = 0; i < dim; ++i)
  {
    double v = values[i];

    if (!sleqp_is_zero(v, zero_eps))
    {
      SLEQP_CALL(sleqp_sparse_vector_push(vec, i, v));
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_sparse_vector_fill(SleqpSparseVec* vec, double value)
{
  if (value == 0)
  {
    SLEQP_CALL(sleqp_sparse_vector_clear(vec));

    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_sparse_vector_reserve(vec, vec->dim));

  SLEQP_CALL(sleqp_sparse_vector_clear(vec));

  for (int k = 0; k < vec->dim; ++k)
  {
    vec->data[k]    = value;
    vec->indices[k] = k;
  }

  vec->nnz = vec->dim;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_sparse_vector_to_raw(const SleqpSparseVec* vec, double* values)
{
  for (int i = 0; i < vec->dim; ++i)
  {
    values[i] = 0.;
  }

  for (int k = 0; k < vec->nnz; ++k)
  {
    values[vec->indices[k]] = vec->data[k];
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_sparse_vector_copy(const SleqpSparseVec* source, SleqpSparseVec* target)
{
  assert(source->dim == target->dim);

  SLEQP_CALL(sleqp_sparse_vector_reserve(target, source->nnz));

  target->nnz = 0;

  for (int k = 0; k < source->nnz; ++k)
  {
    SLEQP_CALL(
      sleqp_sparse_vector_push(target, source->indices[k], source->data[k]));
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_sparse_vector_clear(SleqpSparseVec* vec)
{
  vec->nnz = 0;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_sparse_vector_reserve(SleqpSparseVec* vec, int nnz_max)
{
  if (vec->nnz_max >= nnz_max)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_realloc(&vec->data, nnz_max));
  SLEQP_CALL(sleqp_realloc(&vec->indices, nnz_max));

  vec->nnz_max = nnz_max;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_sparse_vector_resize(SleqpSparseVec* vec, int dim)
{
  if (dim < vec->dim)
  {
    while (vec->nnz > 0 && vec->indices[vec->nnz - 1] >= dim)
    {
      --vec->nnz;
    }
  }

  vec->dim = dim;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_sparse_vector_concat(const SleqpSparseVec* first,
                           const SleqpSparseVec* second,
                           SleqpSparseVec* result)
{
  assert(first != result);
  assert(second != result);

  SLEQP_CALL(sleqp_sparse_vector_clear(result));
  SLEQP_CALL(sleqp_sparse_vector_reserve(result, first->nnz + second->nnz));
  SLEQP_CALL(sleqp_sparse_vector_resize(result, first->dim + second->dim));

  for (int k = 0; k < first->nnz; ++k)
  {
    SLEQP_CALL(
      sleqp_sparse_vector_push(result, first->indices[k], first->data[k]));
  }

  for (int k = 0; k < second->nnz; ++k)
  {
    SLEQP_CALL(sleqp_sparse_vector_push(result,
                                        second->indices[k] + first->dim,
                                        second->data[k]));
  }

  return SLEQP_OKAY;
}

bool
sleqp_sparse_vector_eq(const SleqpSparseVec* first,
                       const SleqpSparseVec* second,
                       double eps)
{
  assert(first->dim == second->dim);

  int k_first = 0, k_second = 0;

  while (k_first < first->nnz || k_second < second->nnz)
  {
    bool valid_first  = (k_first < first->nnz);
    bool valid_second = (k_second < second->nnz);

    int i_first  = valid_first ? first->indices[k_first] : first->dim + 1;
    int i_second = valid_second ? second->indices[k_second] : second->dim + 1;

    int i_combined = SLEQP_MIN(i_first, i_second);

    valid_first  = valid_first && i_first == i_combined;
    valid_second = valid_second && i_second == i_combined;

    double first_value  = valid_first ? first->data[k_first] : 0.;
    double second_value = valid_second ? second->data[k_second] : 0.;

    if (!sleqp_is_eq(first_value, second_value, eps))
    {
      return false;
    }

    if (valid_first)
    {
      ++k_first;
    }

    if (valid_second)
    {
      ++k_second;
    }
  }

  return true;
}

SLEQP_RETCODE
sleqp_sparse_vector_dot(const SleqpSparseVec* first,
                        const SleqpSparseVec* second,
                        double* product)
{
  assert(first->dim == second->dim);

  *product = 0.;

  int k_first = 0, k_second = 0;

  while (k_first < first->nnz && k_second < second->nnz)
  {
    int i_first  = first->indices[k_first];
    int i_second = second->indices[k_second];

    if (i_first == i_second)
    {
      *product += first->data[k_first] * second->data[k_second];
      ++k_first;
      ++k_second;
    }
    else if (i_first < i_second)
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

SLEQP_RETCODE
sleqp_sparse_vector_scale(SleqpSparseVec* vector, const double factor)
{
  // assert(!(isinf(factor) || isnan(factor)));

  if (factor == 0.)
  {
    SLEQP_CALL(sleqp_sparse_vector_clear(vector));
  }
  else if (factor == 1.)
  {
    return SLEQP_OKAY;
  }

  for (int k = 0; k < vector->nnz; ++k)
  {
    vector->data[k] *= factor;
  }

  return SLEQP_OKAY;
}

double
sleqp_sparse_vector_norm(const SleqpSparseVec* vec)
{
  double normsq = sleqp_sparse_vector_norm_sq(vec);

  return sqrt(normsq);
}

double
sleqp_sparse_vector_one_norm(const SleqpSparseVec* vec)
{
  double one_norm = 0.;

  for (int k = 0; k < vec->nnz; ++k)
  {
    double value = SLEQP_ABS(vec->data[k]);

    one_norm += value;
  }

  return one_norm;
}

double
sleqp_sparse_vector_norm_sq(const SleqpSparseVec* vec)
{
  double normsq = 0.;

  for (int k = 0; k < vec->nnz; ++k)
  {
    double value = vec->data[k];

    normsq += value * value;
  }

  return normsq;
}

double
sleqp_sparse_vector_inf_norm(const SleqpSparseVec* vec)
{
  double norminf = 0.;

  for (int k = 0; k < vec->nnz; ++k)
  {
    double value = SLEQP_ABS(vec->data[k]);

    norminf = SLEQP_MAX(value, norminf);
  }

  return norminf;
}

SLEQP_RETCODE
sleqp_sparse_vector_add(const SleqpSparseVec* first,
                        const SleqpSparseVec* second,
                        const double eps,
                        SleqpSparseVec* result)
{
  SLEQP_CALL(
    sleqp_sparse_vector_add_scaled(first, second, 1., 1., eps, result));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_sparse_vector_add_scaled(const SleqpSparseVec* first,
                               const SleqpSparseVec* second,
                               const double first_factor,
                               const double second_factor,
                               const double eps,
                               SleqpSparseVec* result)
{
  assert(first->dim == second->dim);

  // assert(!(isinf(first_factor) || isnan(first_factor)));
  // assert(!(isinf(second_factor) || isnan(second_factor)));

  assert(first != result);
  assert(second != result);

  SLEQP_CALL(sleqp_sparse_vector_clear(result));

  SLEQP_CALL(sleqp_sparse_vector_resize(result, first->dim));

  SLEQP_CALL(sleqp_sparse_vector_reserve(result, first->nnz + second->nnz));

  int k_first = 0, k_second = 0;

  while (k_first < first->nnz || k_second < second->nnz)
  {
    bool valid_first  = (k_first < first->nnz);
    bool valid_second = (k_second < second->nnz);

    int i_first  = valid_first ? first->indices[k_first] : first->dim + 1;
    int i_second = valid_second ? second->indices[k_second] : second->dim + 1;

    int i_combined = SLEQP_MIN(i_first, i_second);

    valid_first  = valid_first && (i_first == i_combined);
    valid_second = valid_second && (i_second == i_combined);

    double value = 0.;

    if (valid_first)
    {
      value += first_factor * first->data[k_first];
      ++k_first;
    }

    if (valid_second)
    {
      value += second_factor * second->data[k_second];
      ++k_second;
    }

    if (sleqp_is_zero(value, eps))
    {
      continue;
    }

    SLEQP_CALL(sleqp_sparse_vector_push(result, i_combined, value));
  }

  return SLEQP_OKAY;
}

double*
sleqp_sparse_vector_at(const SleqpSparseVec* vec, int index)
{
  assert(index < vec->dim);

  for (int i = 0; i < vec->nnz; ++i)
  {
    if (vec->indices[i] == index)
    {
      return vec->data + i;
    }
  }

  return NULL;
}

double
sleqp_sparse_vector_value_at(const SleqpSparseVec* vec, int index)
{
  const double* ptr = sleqp_sparse_vector_at(vec, index);

  return ptr ? (*ptr) : 0.;
}

bool
sleqp_sparse_vector_is_boxed(const SleqpSparseVec* x,
                             const SleqpSparseVec* lb,
                             const SleqpSparseVec* ub)
{
  const int dim = x->dim;

  assert(lb->dim == dim);
  assert(ub->dim == dim);

  int k_x = 0, k_lb = 0, k_ub = 0;

  while (k_x < x->nnz || k_lb < lb->nnz || k_ub < ub->nnz)
  {
    bool valid_x  = (k_x < x->nnz);
    bool valid_lb = (k_lb < lb->nnz);
    bool valid_ub = (k_ub < ub->nnz);

    int idx = valid_x ? x->indices[k_x] : dim + 1;
    idx     = SLEQP_MIN(idx, valid_lb ? lb->indices[k_lb] : dim + 1);
    idx     = SLEQP_MIN(idx, valid_ub ? ub->indices[k_ub] : dim + 1);

    valid_x  = valid_x && idx == x->indices[k_x];
    valid_lb = valid_lb && idx == lb->indices[k_lb];
    valid_ub = valid_ub && idx == ub->indices[k_ub];

    double x_val  = valid_x ? x->data[k_x] : 0.;
    double lb_val = valid_lb ? lb->data[k_lb] : 0.;
    double ub_val = valid_ub ? ub->data[k_ub] : 0.;

    if (x_val < lb_val)
    {
      return false;
    }

    if (x_val > ub_val)
    {
      return false;
    }

    if (valid_lb)
    {
      ++k_lb;
    }

    if (valid_ub)
    {
      ++k_ub;
    }

    if (valid_x)
    {
      ++k_x;
    }
  }

  return true;
}

SLEQP_RETCODE
sleqp_sparse_vector_clip(const SleqpSparseVec* x,
                         const SleqpSparseVec* lb,
                         const SleqpSparseVec* ub,
                         const double eps,
                         SleqpSparseVec* xclip)
{
  const int dim = x->dim;

  assert(x != xclip);
  assert(lb != xclip);
  assert(ub != xclip);

  assert(lb->dim == dim);
  assert(ub->dim == dim);
  assert(xclip->dim == dim);

  SLEQP_CALL(sleqp_sparse_vector_clear(xclip));

  {
    int nnz = SLEQP_MAX(lb->nnz, x->nnz);
    nnz += SLEQP_MAX(ub->nnz, x->nnz);

    nnz = SLEQP_MIN(nnz, dim);

    SLEQP_CALL(sleqp_sparse_vector_reserve(xclip, nnz));
  }

  int k_x = 0, k_lb = 0, k_ub = 0;

  while (k_x < x->nnz || k_lb < lb->nnz || k_ub < ub->nnz)
  {
    bool valid_x  = (k_x < x->nnz);
    bool valid_lb = (k_lb < lb->nnz);
    bool valid_ub = (k_ub < ub->nnz);

    int i_x  = valid_x ? x->indices[k_x] : dim + 1;
    int i_lb = valid_lb ? lb->indices[k_lb] : dim + 1;
    int i_ub = valid_ub ? ub->indices[k_ub] : dim + 1;

    int i_combined;

    i_combined = SLEQP_MIN(i_lb, i_ub);
    i_combined = SLEQP_MIN(i_combined, i_x);

    valid_x  = valid_x && (i_x == i_combined);
    valid_lb = valid_lb && (i_lb == i_combined);
    valid_ub = valid_ub && (i_ub == i_combined);

    double x_val  = valid_x ? x->data[k_x] : 0.;
    double lb_val = valid_lb ? lb->data[k_lb] : 0.;
    double ub_val = valid_ub ? ub->data[k_ub] : 0.;

    x_val = SLEQP_MIN(x_val, ub_val);
    x_val = SLEQP_MAX(x_val, lb_val);

    if (!sleqp_is_zero(x_val, eps))
    {
      SLEQP_CALL(sleqp_sparse_vector_push(xclip, i_combined, x_val));
    }

    if (valid_x)
    {
      ++k_x;
    }

    if (valid_lb)
    {
      ++k_lb;
    }

    if (valid_ub)
    {
      ++k_ub;
    }
  }

  assert(sleqp_sparse_vector_is_boxed(xclip, lb, ub));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_sparse_vector_remove_entries(const SleqpSparseVec* source,
                                   SleqpSparseVec* target,
                                   const int* entry_indices,
                                   int num_entries)
{
  SLEQP_CALL(sleqp_sparse_vector_clear(target));

  assert(source->dim == target->dim + num_entries);

  SLEQP_CALL(sleqp_sparse_vector_reserve(target, source->nnz));

  int k_f    = 0;
  int offset = 0;

  for (int k = 0; k < source->nnz; ++k)
  {
    const int i    = source->indices[k];
    const double v = source->data[k];

    while (k_f < num_entries && entry_indices[k_f] < i)
    {
      ++k_f;
      ++offset;
    }

    if (k_f < num_entries && entry_indices[k_f] == i)
    {
      continue;
    }

    SLEQP_CALL(sleqp_sparse_vector_push(target, i - offset, v));
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_sparse_vector_fprintf(const SleqpSparseVec* vec, FILE* output)
{
  fprintf(output,
          "Sparse vector, dimension: %d, entries: %d\n",
          vec->dim,
          vec->nnz);

  for (int index = 0; index < vec->nnz; ++index)
  {
    fprintf(output, "(%d) = %.14e\n", vec->indices[index], vec->data[index]);
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_sparse_vector_dump(const SleqpSparseVec* vec, FILE* output)
{
  int k = 0;

  for (int i = 0; i < vec->dim; ++i)
  {
    if (k >= vec->nnz)
    {
      fprintf(output, "0.\n");
      continue;
    }

    if (vec->indices[k] == i)
    {
      fprintf(output, "%.14e\n", vec->data[k]);
      ++k;
    }
    else
    {
      fprintf(output, "0.\n");
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_sparse_vector_dump_to_file(const SleqpSparseVec* vec, const char* name)
{
  FILE* output = fopen(name, "w");

  if (!output)
  {
    sleqp_raise(SLEQP_ILLEGAL_ARGUMENT,
                "Failed to open output file '%s'",
                name);
  }

  SLEQP_CALL(sleqp_sparse_vector_dump(vec, output));

  fclose(output);

  return SLEQP_OKAY;
}

bool
sleqp_sparse_vector_is_valid(const SleqpSparseVec* vec)
{
  if (vec->nnz > vec->nnz_max || vec->nnz < 0)
  {
    return false;
  }

  if (vec->nnz == 0)
  {
    return true;
  }

  for (int k = 0; k < vec->nnz; ++k)
  {
    if (vec->indices[k] < 0)
    {
      return false;
    }
  }

  for (int k = 0; k < vec->nnz - 1; ++k)
  {
    if (vec->indices[k] >= vec->indices[k + 1])
    {
      return false;
    }
  }

  if (vec->indices[vec->nnz - 1] >= vec->dim)
  {
    return false;
  }

  return true;
}

bool
sleqp_sparse_vector_is_finite(const SleqpSparseVec* vec)
{
  for (int k = 0; k < vec->nnz - 1; ++k)
  {
    if (isnan(vec->data[k]) || isinf(vec->data[k]))
    {
      return false;
    }
  }

  return true;
}

SLEQP_RETCODE
sleqp_sparse_vector_free(SleqpSparseVec** vstar)
{
  SleqpSparseVec* vec = *vstar;

  if (!vec)
  {
    return SLEQP_OKAY;
  }

  sleqp_free(&vec->indices);
  sleqp_free(&vec->data);

  sleqp_free(&vec);

  *vstar = NULL;

  return SLEQP_OKAY;
}
