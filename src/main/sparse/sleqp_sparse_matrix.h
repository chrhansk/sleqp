#ifndef SLEQP_SPARSE_MATRIX_H
#define SLEQP_SPARSE_MATRIX_H

/**
 * @file sleqp_sparse_matrix.h
 * @brief Definition of sparse matrices.
 **/

#ifdef __cplusplus
extern "C" {
#endif

#include "sleqp_types.h"
#include "sleqp_sparse_vec.h"

  /**
   * A sparse matrix data structure.
   * So far the data is stored in CSC format.
   * Specifically:
   *
   * \f$ A(i, j) = data[k]\f$ iff \f$ rows[k] = i \f$ and \f$ cols[j] <= k < cols[j + 1] \f$
   *
   * for \f$ k = 0, \ldots, nnz - 1 \f$
   *
   **/
  typedef struct SleqpSparseMatrix SleqpSparseMatrix;

  SLEQP_RETCODE sleqp_sparse_matrix_create(SleqpSparseMatrix** matrix,
                                           int num_rows,
                                           int num_cols,
                                           int nnz_max);

  SLEQP_RETCODE sleqp_sparse_matrix_reserve(SleqpSparseMatrix* matrix,
                                            int nnz);

  SLEQP_RETCODE sleqp_sparse_matrix_resize(SleqpSparseMatrix* matrix,
                                           int num_rows,
                                           int num_cols);

  int sleqp_sparse_matrix_get_num_cols(SleqpSparseMatrix* matrix);

  int sleqp_sparse_matrix_get_num_rows(SleqpSparseMatrix* matrix);

  int sleqp_sparse_matrix_get_nnz(SleqpSparseMatrix* matrix);

  int sleqp_sparse_matrix_get_nnz_max(SleqpSparseMatrix* matrix);

  SLEQP_RETCODE sleqp_sparse_matrix_set_nnz(SleqpSparseMatrix* matrix,
                                            int nnz);

  bool sleqp_sparse_matrix_is_quadratic(SleqpSparseMatrix* matrix);


  double* sleqp_sparse_matrix_get_data(SleqpSparseMatrix* matrix);

  int* sleqp_sparse_matrix_get_cols(SleqpSparseMatrix* matrix);

  int* sleqp_sparse_matrix_get_rows(SleqpSparseMatrix* matrix);

  SLEQP_RETCODE sleqp_sparse_matrix_push(SleqpSparseMatrix* matrix,
                                         int row,
                                         int col,
                                         double value);

  SLEQP_RETCODE sleqp_sparse_matrix_push_column(SleqpSparseMatrix* matrix,
                                                int col);

  SLEQP_RETCODE sleqp_sparse_matrix_pop_column(SleqpSparseMatrix* matrix,
                                               int col);

  SLEQP_RETCODE sleqp_sparse_matrix_vector_product(SleqpSparseMatrix* matrix,
                                                   SleqpSparseVec* vector,
                                                   double* result);

  SLEQP_RETCODE sleqp_sparse_matrix_trans_vector_product(SleqpSparseMatrix* matrix,
                                                         SleqpSparseVec* vector,
                                                         double eps,
                                                         SleqpSparseVec* result);

  double* sleqp_sparse_matrix_at(SleqpSparseMatrix* matrix,
                                 int row,
                                 int col);

  bool sleqp_sparse_matrix_eq(SleqpSparseMatrix* first,
                              SleqpSparseMatrix* second,
                              double eps);

  SLEQP_RETCODE sleqp_sparse_matrix_clear(SleqpSparseMatrix* matrix);

  SLEQP_RETCODE sleqp_sparse_matrix_fprintf(SleqpSparseMatrix* matrix,
                                            FILE* output);

  SLEQP_RETCODE sleqp_sparse_matrix_dump(SleqpSparseMatrix* matrix,
                                         FILE* output);

  SLEQP_RETCODE sleqp_sparse_matrix_copy(SleqpSparseMatrix* source,
                                         SleqpSparseMatrix* target);

  bool sleqp_sparse_matrix_valid(SleqpSparseMatrix* matrix);

  SLEQP_RETCODE sleqp_sparse_matrix_capture(SleqpSparseMatrix* matrix);

  SLEQP_RETCODE sleqp_sparse_matrix_release(SleqpSparseMatrix** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_SPARSE_MATRIX_H */
