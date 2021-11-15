#include "preprocessing.h"

#include "fail.h"
#include <assert.h>

SLEQP_RETCODE
sleqp_preprocessing_merge_entries(const SleqpSparseVec* source,
                                  SleqpSparseVec* target,
                                  int num_entries,
                                  const int* entry_indices,
                                  double* entry_values)
{
  SLEQP_CALL(sleqp_sparse_vector_clear(target));

  assert(source->dim + num_entries == target->dim);

  SLEQP_CALL(sleqp_sparse_vector_reserve(target, source->nnz + num_entries));

  int offset = 0;
  int k_f    = 0;

  for (int k = 0; k < source->nnz; ++k)
  {
    const int i_v = source->indices[k];

    while (k_f < num_entries && entry_indices[k_f] <= i_v + offset)
    {
      const int i_f = entry_indices[k_f];

      SLEQP_CALL(sleqp_sparse_vector_push(target, i_f, entry_values[k_f]));

      ++k_f;
      ++offset;
    }

    SLEQP_CALL(sleqp_sparse_vector_push(target, i_v + offset, source->data[k]));
  }

  while (k_f < num_entries)
  {
    SLEQP_CALL(
      sleqp_sparse_vector_push(target, entry_indices[k_f], entry_values[k_f]));

    ++k_f;
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_preprocessing_add_zero_entries(const SleqpSparseVec* source,
                                     SleqpSparseVec* target,
                                     int num_entries,
                                     const int* entry_indices)
{
  SLEQP_CALL(sleqp_sparse_vector_clear(target));

  SLEQP_CALL(sleqp_sparse_vector_reserve(target, source->nnz));

  int k_f    = 0;
  int offset = 0;

  for (int k = 0; k < source->nnz; ++k)
  {
    const int i    = source->indices[k];
    const double v = source->data[k];

    while (k_f < num_entries && entry_indices[k_f] <= i + offset)
    {
      ++k_f;
      ++offset;
    }

    SLEQP_CALL(sleqp_sparse_vector_push(target, i + offset, v));
  }

  return SLEQP_OKAY;
}
