#include "sleqp_preprocessing.h"

SLEQP_RETCODE sleqp_preprocessing_merge_fixed_entries(const SleqpSparseVec* source,
                                                      SleqpSparseVec* target,
                                                      int num_fixed,
                                                      const int* fixed_indices,
                                                      double* fixed_values)
{
  SLEQP_CALL(sleqp_sparse_vector_clear(target));

  int offset = 0;
  int k_f = 0;

  for(int k = 0; k < source->nnz; ++k)
  {
    const int i_v = source->indices[k];

    while(k_f < num_fixed &&
          fixed_indices[k_f] <= i_v)
    {
      const int i_f = fixed_indices[k_f];

      SLEQP_CALL(sleqp_sparse_vector_push(target,
                                          i_f,
                                          fixed_values[k_f]));

      ++k_f;
      ++offset;
    }

    SLEQP_CALL(sleqp_sparse_vector_push(target,
                                        i_v + offset,
                                        source->data[k]));
  }

  return SLEQP_OKAY;
}


SLEQP_RETCODE sleqp_preprocessing_remove_fixed_entries(const SleqpSparseVec* source,
                                                       SleqpSparseVec* target,
                                                       int num_fixed,
                                                       const int* fixed_indices)
{
  SLEQP_CALL(sleqp_sparse_vector_clear(target));

  SLEQP_CALL(sleqp_sparse_vector_reserve(target,
                                         source->nnz));

  int k_f = 0;
  int offset = 0;

  for(int k = 0; k < source->nnz; ++k)
  {
    const int i = source->indices[k];
    const double v = source->data[k];

    while(k_f < num_fixed &&
          fixed_indices[k_f] < i)
    {
      ++k_f;
      ++offset;
    }

    if(k_f < num_fixed &&
       fixed_indices[k_f] == i)
    {
      continue;
    }

    SLEQP_CALL(sleqp_sparse_vector_push(target,
                                        i - offset,
                                        v));
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_preprocessing_remove_fixed_matrix_entries(const SleqpSparseMatrix* source,
                                                              SleqpSparseMatrix* target,
                                                              int num_fixed,
                                                              const int* fixed_indices)
{
  SLEQP_CALL(sleqp_sparse_matrix_reserve(target,
                                         sleqp_sparse_matrix_get_nnz(source)));

  SLEQP_CALL(sleqp_sparse_matrix_clear(target));

  const int num_cols = sleqp_sparse_matrix_get_num_rows(source);

  double* source_data = sleqp_sparse_matrix_get_data(source);
  int* source_rows = sleqp_sparse_matrix_get_rows(source);
  int* source_cols = sleqp_sparse_matrix_get_cols(source);

  int offset = 0;

  for(int col = 0; col < num_cols; ++col)
  {
    SLEQP_CALL(sleqp_sparse_matrix_push_column(target,
                                               col));

    if(offset < num_fixed && fixed_indices[offset] <= col)
    {
      ++offset;
      continue;
    }

    for(int k = source_cols[col]; k < source_cols[col + 1]; ++k)
    {
      SLEQP_CALL(sleqp_sparse_matrix_push(target,
                                          source_rows[k],
                                          col,
                                          source_data[k]));
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_preprocessing_add_zero_entries(const SleqpSparseVec* source,
                                                   SleqpSparseVec* target,
                                                   int num_fixed,
                                                   const int* fixed_indices)
{
  SLEQP_CALL(sleqp_sparse_vector_clear(target));

  SLEQP_CALL(sleqp_sparse_vector_reserve(target,
                                         source->nnz));

  int k_f = 0;
  int offset = 0;

  for(int k = 0; k < source->nnz; ++k)
  {
    const int i = source->indices[k];
    const double v = source->data[k];

    while(k_f < num_fixed && fixed_indices[k] <= i)
    {
      ++k_f;
      ++offset;
    }

    SLEQP_CALL(sleqp_sparse_vector_push(target,
                                        i + offset,
                                        v));
  }

  return SLEQP_OKAY;
}
