#ifndef SLEQP_PUB_SPARSE_MATRIX_H
#define SLEQP_PUB_SPARSE_MATRIX_H

/**
 * @file sparse_matrix.h
 * @brief Definition of sparse matrices.
 **/

#include "sleqp/export.h"
#include "sleqp/pub_types.h"
#include "sleqp/sparse/pub_sparse_vec.h"

/**
 * A sparse matrix data structure.
 * So far the data is stored in CSC format.
 * Specifically:
 *
 * \f$ A(i, j) = data[k]\f$ iff \f$ rows[k] = i \f$ and \f$ cols[j] <= k <
 *cols[j + 1] \f$
 *
 * for \f$ k = 0, \ldots, nnz - 1 \f$
 *
 **/
typedef struct SleqpSparseMatrix SleqpSparseMatrix;

/**
 * Creates a new sparse matrix with a specified number of nonzeros
 *
 * @param[in] matrix     A pointer to the matrix to be created
 * @param[in] num_rows   The desired number of rows
 * @param[in] num_cols   The desired number of columns
 * @param[in] nnz_max    The desired number of nonzeros
 **/
SLEQP_EXPORT SLEQP_NODISCARD SLEQP_RETCODE
sleqp_sparse_matrix_create(SleqpSparseMatrix** matrix,
                           int num_rows,
                           int num_cols,
                           int nnz_max);

/**
 * Reserves a number of nonzeros for the given matrix
 *
 * @param[in] matrix   The matrix
 * @param[in] nnz      The desired number of nonzeros
 **/
SLEQP_EXPORT SLEQP_NODISCARD SLEQP_RETCODE
sleqp_sparse_matrix_reserve(SleqpSparseMatrix* matrix, int nnz);

/**
 * Resizes the given matrix
 *
 * @note If the matrix is non-empty, decreasing the size can leave the matrix in
 *an inconsistent state
 *
 * @param[in] matrix     The matrix
 * @param[in] num_rows   The desired number of rows
 * @param[in] num_cols   The desired number of columns
 **/
SLEQP_EXPORT SLEQP_NODISCARD SLEQP_RETCODE
sleqp_sparse_matrix_resize(SleqpSparseMatrix* matrix,
                           int num_rows,
                           int num_cols);

SLEQP_EXPORT SLEQP_NODISCARD SLEQP_RETCODE
sleqp_sparse_matrix_scale(SleqpSparseMatrix* matrix, double scale);

/**
 * Returns the number of columns of the given matrix
 **/
SLEQP_EXPORT int
sleqp_sparse_matrix_num_cols(const SleqpSparseMatrix* matrix);

/**
 * Returns the number of rows of the given matrix
 **/
SLEQP_EXPORT int
sleqp_sparse_matrix_num_rows(const SleqpSparseMatrix* matrix);

/**
 * Returns the number of nonzeros of the given matrix
 **/
SLEQP_EXPORT int
sleqp_sparse_matrix_nnz(const SleqpSparseMatrix* matrix);

/**
 * Returns the maximum number of nonzeros of the given matrix
 **/
SLEQP_EXPORT int
sleqp_sparse_matrix_nnz_max(const SleqpSparseMatrix* matrix);

/**
 * Sets the number of nonzeros of the given matrix
 **/
SLEQP_EXPORT SLEQP_NODISCARD SLEQP_RETCODE
sleqp_sparse_matrix_set_nnz(SleqpSparseMatrix* matrix, int nnz);

/**
 * Returns whether the given matrix is rectangular
 **/
SLEQP_EXPORT bool
sleqp_sparse_matrix_is_quadratic(const SleqpSparseMatrix* matrix);

/**
 * Returns a pointer to the values of the matrix
 **/
SLEQP_EXPORT double*
sleqp_sparse_matrix_data(const SleqpSparseMatrix* matrix);

/**
 * Returns a pointer to the columns of the given matrix
 **/
SLEQP_EXPORT int*
sleqp_sparse_matrix_cols(const SleqpSparseMatrix* matrix);

/**
 * Returns a pointer to the rows of the given matrix
 **/
SLEQP_EXPORT int*
sleqp_sparse_matrix_rows(const SleqpSparseMatrix* matrix);

/**
 * Pushes a new entry to the matrix. Fails if the matrix is at capacity
 **/
SLEQP_EXPORT SLEQP_NODISCARD SLEQP_RETCODE
sleqp_sparse_matrix_push(SleqpSparseMatrix* matrix,
                         int row,
                         int col,
                         double value);

SLEQP_EXPORT SLEQP_NODISCARD SLEQP_RETCODE
sleqp_sparse_matrix_push_vec(SleqpSparseMatrix* matrix,
                             int col,
                             SleqpSparseVec* vec);

SLEQP_EXPORT SLEQP_NODISCARD SLEQP_RETCODE
sleqp_sparse_matrix_push_column(SleqpSparseMatrix* matrix, int col);

SLEQP_EXPORT SLEQP_NODISCARD SLEQP_RETCODE
sleqp_sparse_matrix_pop_column(SleqpSparseMatrix* matrix, int col);

/**
 * Increases the reference count of the given matrix
 */
SLEQP_EXPORT SLEQP_NODISCARD SLEQP_RETCODE
sleqp_sparse_matrix_capture(SleqpSparseMatrix* matrix);

/**
 * Decreases the reference count of the given matrix, freeing it
 * if the reference count reaches count
 */
SLEQP_EXPORT SLEQP_NODISCARD SLEQP_RETCODE
sleqp_sparse_matrix_release(SleqpSparseMatrix** star);

#endif /* SLEQP_PUB_SPARSE_MATRIX_H */
