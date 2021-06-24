#ifndef SLEQP_PREPROCESSING_H
#define SLEQP_PREPROCESSING_H

#include "sparse/sleqp_sparse_matrix.h"
#include "sparse/sleqp_sparse_vec.h"

#ifdef __cplusplus
extern "C" {
#endif

  SLEQP_RETCODE sleqp_preprocessing_merge_fixed_entries(const SleqpSparseVec* source,
                                                        SleqpSparseVec* target,
                                                        int num_fixed,
                                                        const int* fixed_indices,
                                                        double* fixed_values);

  SLEQP_RETCODE sleqp_preprocessing_remove_fixed_entries(const SleqpSparseVec* source,
                                                         SleqpSparseVec* target,
                                                         int num_fixed,
                                                         const int* fixed_indices);

  SLEQP_RETCODE sleqp_preprocessing_remove_fixed_matrix_entries(const SleqpSparseMatrix* source,
                                                                SleqpSparseMatrix* target,
                                                                int num_fixed,
                                                                const int* fixed_indices);

  SLEQP_RETCODE sleqp_preprocessing_add_zero_entries(const SleqpSparseVec* source,
                                                     SleqpSparseVec* target,
                                                     int num_fixed,
                                                     const int* fixed_indices);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_PREPROCESSING_H */
