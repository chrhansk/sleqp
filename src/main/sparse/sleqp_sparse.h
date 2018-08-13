#ifndef SLEQP_SPARSE_H
#define SLEQP_SPARSE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>

#include "sleqp_types.h"

  /**
   * A sparse vector data structure. Indices
   * are stored in an ascending fashion.
   **/
  typedef struct SleqpSparseVec
  {
    double* data;
    size_t* indices;

    size_t dim;
    size_t nnz;
    size_t nnz_max;

  } SleqpSparseVec;

  SLEQP_RETCODE sleqp_sparse_vector_create(SleqpSparseVec** vec,
                                           size_t dim,
                                           size_t nnz_max);

  SLEQP_RETCODE sleqp_sparse_vector_push(SleqpSparseVec* vec,
                                         size_t idx,
                                         double value);

  SLEQP_RETCODE sleqp_sparse_vector_reserve(SleqpSparseVec* vec,
                                            size_t nnz);

  SLEQP_RETCODE sleqp_sparse_vector_clip(SleqpSparseVec* x,
                                         SleqpSparseVec* lb,
                                         SleqpSparseVec* ub,
                                         SleqpSparseVec** xstar);

  SLEQP_RETCODE sleqp_sparse_vector_free(SleqpSparseVec** vec);

  /**
   * A sparse matrix data structure.
   * So far the data is stored in CSC format.
   **/
  typedef struct SleqpSparseMatrix
  {
    size_t num_rows, num_cols;

    size_t nnz;
    size_t nnz_max;

    double* data;
    size_t* cols;
    size_t* rows;

  } SleqpSparseMatrix;

  SLEQP_RETCODE sleqp_sparse_matrix_create(SleqpSparseMatrix** matrix,
                                           size_t num_rows,
                                           size_t num_cols,
                                           size_t nnz_max);

  SLEQP_RETCODE sleqp_sparse_matrix_reserve(SleqpSparseMatrix* matrix,
                                            size_t nnz);

  SLEQP_RETCODE sleqp_sparse_matrix_fprintf(SleqpSparseMatrix* matrix,
                                            FILE* output);

  SLEQP_RETCODE sleqp_sparse_matrix_free(SleqpSparseMatrix** matrix);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_SPARSE_H */
