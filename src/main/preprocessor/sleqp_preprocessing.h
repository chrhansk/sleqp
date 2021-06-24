#ifndef SLEQP_PREPROCESSING_H
#define SLEQP_PREPROCESSING_H

#include "sparse/sleqp_sparse_matrix.h"
#include "sparse/sleqp_sparse_vec.h"

#ifdef __cplusplus
extern "C" {
#endif

  SLEQP_RETCODE sleqp_preprocessing_merge_entries(const SleqpSparseVec* source,
                                                  SleqpSparseVec* target,
                                                  int num_entries,
                                                  const int* entry_indices,
                                                  double* entry_values);

  SLEQP_RETCODE sleqp_preprocessing_remove_entries(const SleqpSparseVec* source,
                                                   SleqpSparseVec* target,
                                                   int num_entries,
                                                   const int* entry_indices);

  SLEQP_RETCODE sleqp_preprocessing_remove_matrix_cols(const SleqpSparseMatrix* source,
                                                       SleqpSparseMatrix* target,
                                                       int num_entries,
                                                       const int* col_indices);

  SLEQP_RETCODE sleqp_preprocessing_remove_matrix_entries(const SleqpSparseMatrix* source,
                                                          SleqpSparseMatrix* target,
                                                          int num_col_entries,
                                                          const int* col_indices,
                                                          int num_row_entries,
                                                          const int* row_indices);

  SLEQP_RETCODE sleqp_preprocessing_add_zero_entries(const SleqpSparseVec* source,
                                                     SleqpSparseVec* target,
                                                     int num_entries,
                                                     const int* entry_indices);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_PREPROCESSING_H */
