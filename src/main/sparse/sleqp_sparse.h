#ifndef SLEQP_SPARSE_H
#define SLEQP_SPARSE_H

/**
 * @file sleqp_sparse.h
 * @brief Definition of SLEQP sparse vectors / matrices.
 **/

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

  /**
   * Creates a new sparse vector. Data and indices are set to
   * have size of nnz_max.
   *
   * @param[out] vec     A pointer to the vector to be created
   * @param[in]  dim     The desired dimension of the vector
   * @param[in]  nnz_max The desired amount of nonzeros of the vector
   *
   **/
  SLEQP_RETCODE sleqp_sparse_vector_create(SleqpSparseVec** vec,
                                           size_t dim,
                                           size_t nnz_max);

  /**
   * Pushes a new entry on top of a sparse vector. The new
   * entry is assumed to have a larger index than the existing ones.
   *
   * @param[in,out] vec    A pointer to the vector
   * @param[in]     idx    The index of the new entry
   * @param[in]     value  The value of the new entry
   **/
  SLEQP_RETCODE sleqp_sparse_vector_push(SleqpSparseVec* vec,
                                         size_t idx,
                                         double value);

  /**
   * Creates the entries of a sparse vector from a dense
   * vector. The sparse vector will reserve an appropriate
   * number of entries, the dimension will be changed to
   * match that of the dense vector.
   *
   * @param[in,out] vec    A pointer to the vector
   * @param[in]     values A vector of values
   * @param[in]     dim    The dimension of the values input
   **/
  SLEQP_RETCODE sleqp_sparse_vector_from_raw(SleqpSparseVec* vec,
                                             double* values,
                                             size_t dim);

  /**
   * Reserves space for additional nonzeros
   *
   * @param[in,out] vec     A pointer to the vector
   * @param[in]     nnz_max The number of nonzeros
   **/
  SLEQP_RETCODE sleqp_sparse_vector_reserve(SleqpSparseVec* vec,
                                            size_t nnz);

  /**
   * Computes the dot product of two sparse vectors
   *
   * @param[in]  first     A pointer to the first vector
   * @param[in]  second    A pointer to the second vector
   * @param[out] product   A pointer to the result
   **/
  SLEQP_RETCODE sleqp_sparse_vector_dot(SleqpSparseVec* first,
                                        SleqpSparseVec* second,
                                        double* product);

  /**
   * Scales the sparse vector by a factor
   *
   * @param[in,out] vector   A pointer to the vector
   * @param[out]    factor   The factor
   **/
  SLEQP_RETCODE sleqp_sparse_vector_scale(SleqpSparseVec* vector,
                                          double factor);

  /**
   * Computes the sum of two sparse vectors
   *
   * @param[in]  first         A pointer to the first vector
   * @param[in]  second        A pointer to the second vector
   * @param[in]  first_factor  A factor for the first vector
   * @param[in]  second_factor A factor for the first vector
   * @param[out] result        A pointer to the result
   **/
  SLEQP_RETCODE sleqp_sparse_vector_add(SleqpSparseVec* first,
                                        SleqpSparseVec* second,
                                        double first_factor,
                                        double second_factor,
                                        SleqpSparseVec* result);

  /**
   * Computes the dot product of a sparse and a dense vector
   *
   * @param[in]  first     A pointer to the first vector
   * @param[in]  second    A pointer to the second vector
   * @param[out] product   A pointer to the result
   **/
  SLEQP_RETCODE sleqp_sparse_vector_dense_dot(SleqpSparseVec* first,
                                              double* second,
                                              double* product);

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
   * \f$ A(i, j) = data[k]\f$ iff \f$ row[k] = i \f$ and \f$ col[j] <= k < col[j + 1] \f$
   *
   * for \f$ k = 0, \ldots, nnz - 1 \f$
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

  SLEQP_RETCODE sleqp_sparse_matrix_vector_product(SleqpSparseMatrix* matrix,
                                                   SleqpSparseVec* vector,
                                                   double* result);

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
