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

  /**
   * Creates a new sparse matrix with a specified number of nonzeros
   *
   * @param[in] matrix     A pointer to the vector to be created
   * @param[in] num_rows   The desired number of rows
   * @param[in] num_cols   The desired number of columns
   * @param[in] nnz_max    The desired number of nonzeros
   **/
  SLEQP_EXPORT SLEQP_RETCODE sleqp_sparse_matrix_create(SleqpSparseMatrix** matrix,
                                                        int num_rows,
                                                        int num_cols,
                                                        int nnz_max);

  /**
   * Reserves a number of nonzeros for the given matrix
   *
   * @param[in] matrix   The matrix
   * @param[in] nnz      The desired number of nonzeros
   **/
  SLEQP_EXPORT SLEQP_RETCODE sleqp_sparse_matrix_reserve(SleqpSparseMatrix* matrix,
                                                         int nnz);

  /**
   * Resizes the given matrix
   *
   * @note If the matrix is non-empty, decreasing the size can leave the matrix in an inconsistent state
   *
   * @param[in] matrix     The matrix
   * @param[in] num_rows   The desired number of rows
   * @param[in] num_cols   The desired number of columns
   **/
  SLEQP_EXPORT SLEQP_RETCODE sleqp_sparse_matrix_resize(SleqpSparseMatrix* matrix,
                                                        int num_rows,
                                                        int num_cols);

  /**
   * Returns the number of columns of the given matrix
   **/
  SLEQP_EXPORT int sleqp_sparse_matrix_get_num_cols(const SleqpSparseMatrix* matrix);

  /**
   * Returns the number of rows of the given matrix
   **/
  SLEQP_EXPORT int sleqp_sparse_matrix_get_num_rows(const SleqpSparseMatrix* matrix);

  /**
   * Returns the number of nonzeros of the given matrix
   **/
  SLEQP_EXPORT int sleqp_sparse_matrix_get_nnz(const SleqpSparseMatrix* matrix);

  /**
   * Returns the maximum number of nonzeros of the given matrix
   **/
  SLEQP_EXPORT int sleqp_sparse_matrix_get_nnz_max(const SleqpSparseMatrix* matrix);

  /**
   * Sets the number of nonzeros of the given matrix
   **/
  SLEQP_EXPORT SLEQP_RETCODE sleqp_sparse_matrix_set_nnz(SleqpSparseMatrix* matrix,
                                                         int nnz);

  /**
   * Returns whether the given matrix is rectangular
   **/
  SLEQP_EXPORT bool sleqp_sparse_matrix_is_quadratic(const SleqpSparseMatrix* matrix);

  /**
   * Returns a pointer to the values of the matrix
   **/
  SLEQP_EXPORT double* sleqp_sparse_matrix_get_data(SleqpSparseMatrix* matrix);

  /**
   * Returns a pointer to the columns of the given matrix
   **/
  SLEQP_EXPORT int* sleqp_sparse_matrix_get_cols(SleqpSparseMatrix* matrix);

  /**
   * Returns a pointer to the rows of the given matrix
   **/
  SLEQP_EXPORT int* sleqp_sparse_matrix_get_rows(SleqpSparseMatrix* matrix);

  /**
   * Pushes a new entry to the matrix. Fails if the matrix is at capacity
   **/
  SLEQP_EXPORT SLEQP_RETCODE sleqp_sparse_matrix_push(SleqpSparseMatrix* matrix,
                                                      int row,
                                                      int col,
                                                      double value);

  SLEQP_EXPORT SLEQP_RETCODE sleqp_sparse_matrix_push_column(SleqpSparseMatrix* matrix,
                                                             int col);

  SLEQP_EXPORT SLEQP_RETCODE sleqp_sparse_matrix_pop_column(SleqpSparseMatrix* matrix,
                                                            int col);

  /**
   * Computes the product of the given matrix with the given vector
   * @param[in]  matrix     The matrix
   * @param[in]  vector     The input vector (dimension equal to the number of rows of the matrix)
   * @param[out] result     The result array (size equal to the number of columns of the matrix)
   **/
  SLEQP_EXPORT SLEQP_RETCODE sleqp_sparse_matrix_vector_product(const SleqpSparseMatrix* matrix,
                                                                const SleqpSparseVec* vector,
                                                                double* result);

  /**
   * Computes the product of the transposed matrix with the given vector
   *
   * @param[in]  matrix     The matrix
   * @param[in]  vector     The input vector (dimension equal to the number of columns of the matrix)
   * @param[in]  eps        The given tolerance
   * @param[out] result     The result vector (dimension equal to the number of rows of the matrix)
   **/
  SLEQP_EXPORT SLEQP_RETCODE sleqp_sparse_matrix_trans_vector_product(const SleqpSparseMatrix* matrix,
                                                                      const SleqpSparseVec* vector,
                                                                      double eps,
                                                                      SleqpSparseVec* result);

  /**
   * Returns a pointer to the entry to the given element of the matrix
   *
   * @param[in] matrix     The matrix
   * @param[in] row        The given row
   * @param[in] col        The given column
   *
   * @return A pointer to the entry, or NULL if the entry is not contained
   **/
  SLEQP_EXPORT double* sleqp_sparse_matrix_at(SleqpSparseMatrix* matrix,
                                              int row,
                                              int col);

  SLEQP_RETCODE sleqp_sparse_lower_triangular(const SleqpSparseMatrix* source,
                                              SleqpSparseMatrix* target);

  /**
   * Returns whether all entries of the given matrices are equal to
   * within the specified tolerance
   *
   * @param[in] first      The first matrix
   * @param[in] second     The second matrix
   * @param[in] eps        The tolerance
   *
   * @sa sleqp_is_eq(double x, double y, double eps)
   **/
  SLEQP_EXPORT bool sleqp_sparse_matrix_eq(const SleqpSparseMatrix* first,
                                           const SleqpSparseMatrix* second,
                                           double eps);

  /**
   * Clears the given matrix, i.e., removes all of its elements
   *
   * @param[in] matrix     The matrix
   **/
  SLEQP_EXPORT SLEQP_RETCODE sleqp_sparse_matrix_clear(SleqpSparseMatrix* matrix);

  /**
   * Prints the matrix to the given file
   *
   * @param[in] matrix     The matrix
   * @param[in] output     The output file
   **/
  SLEQP_EXPORT SLEQP_RETCODE sleqp_sparse_matrix_fprintf(const SleqpSparseMatrix* matrix,
                                                         FILE* output);

  /**
   * Dumps the matrix to the given file using the MatrixMarket format
   *
   * @param[in] matrix     The matrix
   * @param[in] output     The output file
   **/
  SLEQP_EXPORT SLEQP_RETCODE sleqp_sparse_matrix_dump(const SleqpSparseMatrix* matrix,
                                                      FILE* output);

  /**
   * Copies the given matrix
   *
   * @param[in]     source     The source matrix
   * @param[in,out] target     The target matrix
   **/
  SLEQP_EXPORT SLEQP_RETCODE sleqp_sparse_matrix_copy(const SleqpSparseMatrix* source,
                                                      SleqpSparseMatrix* target);

  /**
   * Returns whether the given matrix is valid, i.e. whether the
   * number of rows / columns is non-negative, all rows / columns are
   * within their respective bounds and properly ordered and all
   * entries are free of infs / nans.
   **/
  SLEQP_EXPORT bool sleqp_sparse_matrix_valid(const SleqpSparseMatrix* matrix);

  /**
   * Increases the reference count of the given matrix
   */
  SLEQP_EXPORT SLEQP_RETCODE sleqp_sparse_matrix_capture(SleqpSparseMatrix* matrix);

  /**
   * Decreases the reference count of the given matrix, freeing it
   * if the reference count reaches count
   */
  SLEQP_EXPORT SLEQP_RETCODE sleqp_sparse_matrix_release(SleqpSparseMatrix** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_SPARSE_MATRIX_H */
