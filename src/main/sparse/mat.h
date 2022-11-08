#ifndef SLEQP_MAT_H
#define SLEQP_MAT_H

#include "pub_mat.h"

#include "vec.h"

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_mat_vstack(const SleqpMat* first,
                 const SleqpMat* second,
                 SleqpMat* result);

/**
 * Computes the product of the given matrix with the given vector
 * @param[in]  matrix     The matrix
 * @param[in]  vector     The input vector (dimension equal to the number of
 *rows of the matrix)
 * @param[out] result     The result array (size equal to the number of columns
 *of the matrix)
 **/
SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_mat_mult_vec(const SleqpMat* matrix,
                   const SleqpVec* vector,
                   double* result);

/**
 * Computes the product of the transposed matrix with the given vector
 *
 * @param[in]  matrix     The matrix
 * @param[in]  vector     The input vector (dimension equal to the number of
 *columns of the matrix)
 * @param[in]  eps        The given tolerance
 * @param[out] result     The result vector (dimension equal to the number of
 *rows of the matrix)
 **/
SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_mat_mult_vec_trans(const SleqpMat* matrix,
                         const SleqpVec* vector,
                         double eps,
                         SleqpVec* result);

/**
 * Returns a pointer to the entry to the given element of the matrix
 *
 * @param[in] matrix     The matrix
 * @param[in] row        The given row
 * @param[in] col        The given column
 *
 * @return A pointer to the entry, or NULL if the entry is not contained
 **/
double*
sleqp_mat_at(SleqpMat* matrix, int row, int col);

/**
 * Returns the value of the entry at the given element of the matrix
 *
 * @param[in] matrix     The matrix
 * @param[in] row        The given row
 * @param[in] col        The given column
 *
 * @return A pointer to the entry, or NULL if the entry is not contained
 **/
double
sleqp_mat_value_at(SleqpMat* matrix, int row, int col);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_mat_col(const SleqpMat* matrix, int col, SleqpVec* vec);

bool
sleqp_mat_is_lower(const SleqpMat* matrix);

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
bool
sleqp_mat_eq(const SleqpMat* first, const SleqpMat* second, double eps);

SLEQP_RETCODE
sleqp_mat_remove_rows(const SleqpMat* source,
                      SleqpMat* target,
                      const int* row_indices,
                      int num_row_entries);

SLEQP_RETCODE
sleqp_mat_remove_cols(const SleqpMat* source,
                      SleqpMat* target,
                      const int* col_indices,
                      int num_col_entries);

SLEQP_RETCODE
sleqp_mat_remove_entries(const SleqpMat* source,
                         SleqpMat* target,
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
SLEQP_RETCODE
sleqp_mat_clear(SleqpMat* matrix);

/**
 * Prints the matrix to the given file
 *
 * @param[in] matrix     The matrix
 * @param[in] output     The output file
 **/
SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_mat_fprintf(const SleqpMat* matrix, FILE* output);

/**
 * Dumps the matrix to the given file using the MatrixMarket format
 *
 * @param[in] matrix     The matrix
 * @param[in] output     The output file
 **/
SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_mat_dump(const SleqpMat* matrix, FILE* output);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_mat_dump_to_file(const SleqpMat* matrix, const char* name);

/**
 * Returns whether the given matrix is valid, i.e. whether the
 * number of rows / columns is non-negative, all rows / columns are
 * within their respective bounds and properly ordered.
 **/
bool
sleqp_mat_is_valid(const SleqpMat* matrix);

/**
 * Transposes the given matrix
 *
 * @param[in]  source    The source matrix
 * @param[out] target    The target to store the transposed source
 * @param[in]  row_cache Cache having at least as many entries as rows of source
 **/
SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_mat_trans(const SleqpMat* source, SleqpMat* target, int* row_cache);

/**
 * Returns whether the entries of the given matrix are finite with respect to
 *  \ref sleqp_is_finite(double)
 **/
bool
sleqp_mat_is_finite(const SleqpMat* matrix);

#endif /* SLEQP_MAT_H */
