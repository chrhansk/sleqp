#ifndef SLEQP_SPARSE_VEC_H
#define SLEQP_SPARSE_VEC_H

#ifdef __cplusplus
extern "C" {
#endif

#include "sleqp_types.h"

  /**
   * A sparse vector data structure. Indices
   * are stored in an ascending fashion.
   **/
  typedef struct SleqpSparseVec
  {
    double* data;
    int* indices;

    int dim;
    int nnz;
    int nnz_max;

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
                                           int dim,
                                           int nnz_max);

  /**
   * Pushes a new entry on top of a sparse vector. The new
   * entry is assumed to have a larger index than the existing ones.
   *
   * @param[in,out] vec    A pointer to the vector
   * @param[in]     idx    The index of the new entry
   * @param[in]     value  The value of the new entry
   **/
  SLEQP_RETCODE sleqp_sparse_vector_push(SleqpSparseVec* vec,
                                         int idx,
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
                                             int dim);

  SLEQP_RETCODE sleqp_sparse_vector_to_raw(SleqpSparseVec* vec,
                                           double* values);

  SLEQP_RETCODE sleqp_sparse_vector_copy(SleqpSparseVec* source,
                                         SleqpSparseVec* target);

  /**
   * Reserves space for additional nonzeros
   *
   * @param[in,out] vec     A pointer to the vector
   * @param[in]     nnz_max The number of nonzeros
   **/
  SLEQP_RETCODE sleqp_sparse_vector_reserve(SleqpSparseVec* vec,
                                            int nnz);

  SLEQP_RETCODE sleqp_sparse_vector_resize(SleqpSparseVec* vec,
                                           int dim);

  SLEQP_Bool sleqp_sparse_vector_eq(SleqpSparseVec* first,
                                    SleqpSparseVec* second);

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

  double sleqp_sparse_vector_normsq(SleqpSparseVec* vec);

  double* sleqp_sparse_vector_at(SleqpSparseVec* vec,
                                 int index);

  SLEQP_RETCODE sleqp_sparse_vector_clip(SleqpSparseVec* x,
                                         SleqpSparseVec* lb,
                                         SleqpSparseVec* ub,
                                         SleqpSparseVec** xstar);

  SLEQP_RETCODE sleqp_sparse_vector_fprintf(SleqpSparseVec* vec,
                                            FILE* output);

  SLEQP_Bool sleqp_sparse_vector_valid(SleqpSparseVec* vec);

  SLEQP_RETCODE sleqp_sparse_vector_free(SleqpSparseVec** vec);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_SPARSE_VEC_H */
