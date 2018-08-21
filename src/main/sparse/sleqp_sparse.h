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

  SLEQP_RETCODE sleqp_sparse_vector_from_raw(SleqpSparseVec* vec,
                                             double* values,
                                             size_t dim);

  SLEQP_RETCODE sleqp_sparse_vector_reserve(SleqpSparseVec* vec,
                                            size_t nnz);

  double* sleqp_sparse_vector_at(SleqpSparseVec* vec,
                                 size_t index);

  SLEQP_RETCODE sleqp_sparse_vector_clip(SleqpSparseVec* x,
                                         SleqpSparseVec* lb,
                                         SleqpSparseVec* ub,
                                         SleqpSparseVec** xstar);

  SLEQP_RETCODE sleqp_sparse_vector_fprintf(SleqpSparseVec* vec,
                                            FILE* output);

  SLEQP_RETCODE sleqp_sparse_vector_free(SleqpSparseVec** vec);

  /**
   * A sparse matrix data structure.
   * So far the data is stored in CSC format.
   * Specifically:
   *
   * A(i, j) = data[k] iff row[k] = i and col[j] <= k < col[j + 1]
   *
   * for k = 0, ..., nnz - 1
   *
   **/
  typedef struct SleqpSparseMatrix
  {
    int num_rows, num_cols;

    int nnz;
    int nnz_max;

    double* data;
    int* cols;
    int* rows;

  } SleqpSparseMatrix;

  SLEQP_RETCODE sleqp_sparse_matrix_create(SleqpSparseMatrix** matrix,
                                           size_t num_rows,
                                           size_t num_cols,
                                           size_t nnz_max);

  SLEQP_RETCODE sleqp_sparse_matrix_reserve(SleqpSparseMatrix* matrix,
                                            size_t nnz);

  SLEQP_RETCODE sleqp_sparse_matrix_resize(SleqpSparseMatrix* matrix,
                                           size_t num_rows,
                                           size_t num_cols);

  SLEQP_RETCODE sleqp_sparse_matrix_push(SleqpSparseMatrix* matrix,
                                         size_t row,
                                         size_t col,
                                         double value);

  SLEQP_RETCODE sleqp_sparse_matrix_add_column(SleqpSparseMatrix* matrix,
                                               size_t col);

  SLEQP_RETCODE sleqp_sparse_matrix_remove_column(SleqpSparseMatrix* matrix,
                                                  size_t col);

  double* sleqp_sparse_matrix_at(SleqpSparseMatrix* matrix,
                                 size_t row,
                                 size_t col);

  SLEQP_RETCODE sleqp_sparse_matrix_fprintf(SleqpSparseMatrix* matrix,
                                            FILE* output);

  SLEQP_RETCODE sleqp_sparse_matrix_free(SleqpSparseMatrix** matrix);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_SPARSE_H */
