#ifndef SLEQP_SPARSE_MATRIX_H
#define SLEQP_SPARSE_MATRIX_H

#include "pub_sparse_matrix.h"

#include "sparse_vec.h"

#ifdef __cplusplus
extern "C" {
#endif

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_sparse_matrix_vstack(const SleqpSparseMatrix* first,
                                           const SleqpSparseMatrix* second,
                                           SleqpSparseMatrix* result);

  /**
   * Computes the product of the given matrix with the given vector
   * @param[in]  matrix     The matrix
   * @param[in]  vector     The input vector (dimension equal to the number of rows of the matrix)
   * @param[out] result     The result array (size equal to the number of columns of the matrix)
   **/
  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_sparse_matrix_vector_product(const SleqpSparseMatrix* matrix,
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
  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_sparse_matrix_trans_vector_product(const SleqpSparseMatrix* matrix,
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
  double* sleqp_sparse_matrix_at(SleqpSparseMatrix* matrix,
                                 int row,
                                 int col);

  /**
   * Returns the value of the entry at the given element of the matrix
   *
   * @param[in] matrix     The matrix
   * @param[in] row        The given row
   * @param[in] col        The given column
   *
   * @return A pointer to the entry, or NULL if the entry is not contained
   **/
  double sleqp_sparse_matrix_value_at(SleqpSparseMatrix* matrix,
                                      int row,
                                      int col);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_sparse_matrix_get_col(const SleqpSparseMatrix* matrix,
                                            int col,
                                            SleqpSparseVec* vec);

  SLEQP_NODISCARD
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
  bool sleqp_sparse_matrix_eq(const SleqpSparseMatrix* first,
                              const SleqpSparseMatrix* second,
                              double eps);

  SLEQP_RETCODE sleqp_sparse_matrix_remove_rows(const SleqpSparseMatrix* source,
                                                SleqpSparseMatrix* target,
                                                const int* row_indices,
                                                int num_row_entries);

  SLEQP_RETCODE sleqp_sparse_matrix_remove_cols(const SleqpSparseMatrix* source,
                                                SleqpSparseMatrix* target,
                                                const int* col_indices,
                                                int num_col_entries);

  SLEQP_RETCODE sleqp_sparse_matrix_remove_entries(const SleqpSparseMatrix* source,
                                                   SleqpSparseMatrix* target,
                                                   const int* col_indices,
                                                   int num_col_entries,
                                                   const int* row_indices,
                                                   int num_row_entries);

  /**
   * Clears the given matrix, i.e., removes all of its elements
   *
   * @param[in] matrix     The matrix
   **/
  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_sparse_matrix_clear(SleqpSparseMatrix* matrix);

  /**
   * Prints the matrix to the given file
   *
   * @param[in] matrix     The matrix
   * @param[in] output     The output file
   **/
  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_sparse_matrix_fprintf(const SleqpSparseMatrix* matrix,
                                            FILE* output);

  /**
   * Dumps the matrix to the given file using the MatrixMarket format
   *
   * @param[in] matrix     The matrix
   * @param[in] output     The output file
   **/
  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_sparse_matrix_dump(const SleqpSparseMatrix* matrix,
                                         FILE* output);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_sparse_matrix_dump_to_file(const SleqpSparseMatrix* matrix,
                                                 const char* name);

  /**
   * Copies the given matrix
   *
   * @param[in]     source     The source matrix
   * @param[in,out] target     The target matrix
   **/
  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_sparse_matrix_copy(const SleqpSparseMatrix* source,
                                         SleqpSparseMatrix* target);

  /**
   * Returns whether the given matrix is valid, i.e. whether the
   * number of rows / columns is non-negative, all rows / columns are
   * within their respective bounds and properly ordered.
   **/
  bool sleqp_sparse_matrix_is_valid(const SleqpSparseMatrix* matrix);

  /**
   * Returns whether the entries of the given matrix are finite with respect to
   *  \ref sleqp_is_finite(double)
   **/
  bool sleqp_sparse_matrix_is_finite(const SleqpSparseMatrix* matrix);


#ifdef __cplusplus
}
#endif

#endif /* SLEQP_SPARSE_MATRIX_H */
