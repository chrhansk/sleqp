#include "sleqp_preprocessing.h"

#include "sleqp_assert.h"

SLEQP_RETCODE sleqp_preprocessing_merge_entries(const SleqpSparseVec* source,
                                                      SleqpSparseVec* target,
                                                      int num_entries,
                                                      const int* entry_indices,
                                                      double* entry_values)
{
  SLEQP_CALL(sleqp_sparse_vector_clear(target));

  int offset = 0;
  int k_f = 0;

  for(int k = 0; k < source->nnz; ++k)
  {
    const int i_v = source->indices[k];

    while(k_f < num_entries &&
          entry_indices[k_f] <= i_v)
    {
      const int i_f = entry_indices[k_f];

      SLEQP_CALL(sleqp_sparse_vector_push(target,
                                          i_f,
                                          entry_values[k_f]));

      ++k_f;
      ++offset;
    }

    SLEQP_CALL(sleqp_sparse_vector_push(target,
                                        i_v + offset,
                                        source->data[k]));
  }

  return SLEQP_OKAY;
}


SLEQP_RETCODE sleqp_preprocessing_remove_entries(const SleqpSparseVec* source,
                                                 SleqpSparseVec* target,
                                                 int num_entries,
                                                 const int* entry_indices)
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

    while(k_f < num_entries &&
          entry_indices[k_f] < i)
    {
      ++k_f;
      ++offset;
    }

    if(k_f < num_entries &&
       entry_indices[k_f] == i)
    {
      continue;
    }

    SLEQP_CALL(sleqp_sparse_vector_push(target,
                                        i - offset,
                                        v));
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_preprocessing_remove_matrix_cols(const SleqpSparseMatrix* source,
                                                     SleqpSparseMatrix* target,
                                                     int num_entries,
                                                     const int* col_indices)
{
  SLEQP_CALL(sleqp_preprocessing_remove_matrix_entries(source,
                                                       target,
                                                       num_entries,
                                                       col_indices,
                                                       0,
                                                       NULL));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_preprocessing_remove_matrix_entries(const SleqpSparseMatrix* source,
                                                        SleqpSparseMatrix* target,
                                                        int num_col_entries,
                                                        const int* col_indices,
                                                        int num_row_entries,
                                                        const int* row_indices)
{
  SLEQP_CALL(sleqp_sparse_matrix_reserve(target,
                                         sleqp_sparse_matrix_get_nnz(source)));

  SLEQP_CALL(sleqp_sparse_matrix_clear(target));

  sleqp_assert(sleqp_sparse_matrix_get_num_cols(source) ==
               sleqp_sparse_matrix_get_num_cols(target) + num_col_entries);

  sleqp_assert(sleqp_sparse_matrix_get_num_rows(source) ==
               sleqp_sparse_matrix_get_num_rows(target) + num_row_entries);


  const int num_cols = sleqp_sparse_matrix_get_num_cols(source);

  double* source_data = sleqp_sparse_matrix_get_data(source);
  int* source_rows = sleqp_sparse_matrix_get_rows(source);
  int* source_cols = sleqp_sparse_matrix_get_cols(source);

  int col_offset = 0;

  for(int col = 0; col < num_cols; ++col)
  {
    if(col_offset < num_col_entries && col_indices[col_offset] <= col)
    {
      ++col_offset;
      continue;
    }

    SLEQP_CALL(sleqp_sparse_matrix_push_column(target,
                                               col - col_offset));

    int row_offset = 0;

    for(int k = source_cols[col]; k < source_cols[col + 1]; ++k)
    {
      const int row = source_rows[k];

      while(row_offset < num_row_entries && row_indices[row_offset] < row)
      {
        ++row_offset;
      }

      if(row_offset < num_row_entries && row_indices[row_offset] == row)
      {
        continue;
      }

      SLEQP_CALL(sleqp_sparse_matrix_push(target,
                                          row - row_offset,
                                          col - col_offset,
                                          source_data[k]));
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_preprocessing_add_zero_entries(const SleqpSparseVec* source,
                                                   SleqpSparseVec* target,
                                                   int num_entries,
                                                   const int* entry_indices)
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

    while(k_f < num_entries && entry_indices[k] <= i)
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
