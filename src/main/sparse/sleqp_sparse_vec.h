#ifndef SLEQP_SPARSE_VEC_H
#define SLEQP_SPARSE_VEC_H

/**
 * @file sleqp_sparse_vec.h
 * @brief Definition of sparse vectors.
 **/

#include "sleqp_export.h"
#include "sleqp_types.h"

#ifdef __cplusplus
extern "C" {
#endif

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
   * @param[in]  vec     A pointer to the vector to be created
   * @param[in]  dim     The desired dimension of the vector
   * @param[in]  nnz_max The desired amount of nonzeros of the vector
   *
   **/
  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_sparse_vector_create(SleqpSparseVec** vec,
                                           int dim,
                                           int nnz_max);

  /**
   * Creates a new sparse vector without allocating memory
   * for non-zero entries
   *
   * @param[in]  vec     A pointer to the vector to be created
   * @param[in]  dim     The desired dimension of the vector
   *
   **/
  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_sparse_vector_create_empty(SleqpSparseVec** vec,
                                                 int dim);

  /**
   * Creates a new sparse vector, allocating memory sufficient
   * for `dim` non-zero entries
   *
   * @param[in]  vec     A pointer to the vector to be created
   * @param[in]  dim     The desired dimension of the vector
   *
   **/
  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_sparse_vector_create_full(SleqpSparseVec** vec,
                                                int dim);

  /**
   * Pushes a new entry on top of a sparse vector. The new
   * entry is assumed to have a larger index than the existing ones.
   *
   * @param[in,out] vec    A pointer to the vector
   * @param[in]     idx    The index of the new entry
   * @param[in]     value  The value of the new entry
   **/
  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_sparse_vector_push(SleqpSparseVec* vec,
                                         int idx,
                                         double value);

  /**
   * Creates the entries of a sparse vector from a dense
   * vector. The sparse vector will reserve an appropriate
   * number of entries, the dimension will be changed to
   * match that of the dense vector.
   *
   * @param[in,out] vec         A pointer to the vector
   * @param[in]     values      A vector of values
   * @param[in]     dim         The dimension of the values input
   * @param[in]     zero_eps    The numerical tolerance
   **/
  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_sparse_vector_from_raw(SleqpSparseVec* vec,
                                             double* values,
                                             int dim,
                                             double zero_eps);

  /**
   * Writes the content of this vector into an array. The
   * array is assumed to have a size of at least the dimension
   * of the given vector
   *
   * @param[in] vec     A pointer to the vector
   * @param[in] values  A pointer to the output array
   *
   **/
  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_sparse_vector_to_raw(const SleqpSparseVec* vec,
                                           double* values);

  /**
   * Copies one vector to another
   *
   * @param[in] source     A pointer to the copy source
   * @param[out] target    A pointer to the copy target
   **/
  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_sparse_vector_copy(const SleqpSparseVec* source,
                                         SleqpSparseVec* target);

  /**
   * Clears the given vector, discarding all entries while
   * keeping the dimension constant
   *
   * @param[in,out] vec     A pointer to the vector
   **/
  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_sparse_vector_clear(SleqpSparseVec* vec);

  /**
   * Reserves space for additional nonzeros
   *
   * @param[in,out] vec     A pointer to the vector
   * @param[in]     nnz_max The number of nonzeros
   **/
  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_sparse_vector_reserve(SleqpSparseVec* vec,
                                            int nnz);

  /**
   * Resizes the vector to the given dimension, discarding
   * entries if necessary
   *
   * @param[in,out] vec     A pointer to the vector
   * @param[in]     dim     The new dimension
   **/
  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_sparse_vector_resize(SleqpSparseVec* vec,
                                           int dim);

  /**
   * Returns whether all entries of the given vector are equal
   * up to the given tolerance
   *
   * @param[in]  first     A pointer to the first vector
   * @param[in]  second    A pointer to the second vector
   * @param[in]  eps       The desred tolerance
   *
   * @sa sleqp_is_eq(double x, double y, double eps)
   **/
  SLEQP_EXPORT bool sleqp_sparse_vector_eq(const SleqpSparseVec* first,
                                           const SleqpSparseVec* second,
                                           double eps);

  /**
   * Computes the dot product of two sparse vectors
   *
   * @param[in]  first     A pointer to the first vector
   * @param[in]  second    A pointer to the second vector
   * @param[out] product   A pointer to the result
   **/
  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_sparse_vector_dot(const SleqpSparseVec* first,
                                        const SleqpSparseVec* second,
                                        double* product);

  /**
   * Scales the sparse vector by a factor
   *
   * @param[in,out] vector   A pointer to the vector
   * @param[in]     factor   The factor
   **/
  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_sparse_vector_scale(SleqpSparseVec* vector,
                                          const double factor);

  /**
   * Computes the sum of two sparse vectors
   *
   * @param[in]  first         A pointer to the first vector
   * @param[in]  second        A pointer to the second vector
   * @param[out] result        A pointer to the result
   **/
  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_sparse_vector_add(const SleqpSparseVec* first,
                                        const SleqpSparseVec* second,
                                        const double eps,
                                        SleqpSparseVec* result);

  /**
   * Computes the weighted sum of two sparse vectors
   *
   * @param[in]  first         A pointer to the first vector
   * @param[in]  second        A pointer to the second vector
   * @param[in]  first_factor  A factor for the first vector
   * @param[in]  second_factor A factor for the first vector
   * @param[out] result        A pointer to the result
   **/
  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_sparse_vector_add_scaled(const SleqpSparseVec* first,
                                               const SleqpSparseVec* second,
                                               const double first_factor,
                                               const double second_factor,
                                               const double eps,
                                               SleqpSparseVec* result);

  /**
   * Returns the 2-norm of the given vector
   *
   * @param[in] vector   A pointer to the vector
   **/
  SLEQP_EXPORT double sleqp_sparse_vector_norm(const SleqpSparseVec* vec);

  /**
   * Returns the 1-norm of the given vector
   *
   * @param[in] vector   A pointer to the vector
   **/
  SLEQP_EXPORT double sleqp_sparse_vector_one_norm(const SleqpSparseVec* vec);

  /**
   * Returns the squared 2-norm of the given vector
   *
   * @param[in] vector   A pointer to the vector
   **/
  SLEQP_EXPORT double sleqp_sparse_vector_norm_sq(const SleqpSparseVec* vec);

  /**
   * Returns the oo-norm of the given vector
   *
   * @param[in] vector   A pointer to the vector
   **/
  SLEQP_EXPORT double sleqp_sparse_vector_inf_norm(const SleqpSparseVec* vec);

  /**
   * Returns a pointer to the entry of the given vector at
   * the given index, or `NULL` if the entry is not present.
   *
   * @param[in] vector   A pointer to the vector
   * @param[in] index    The desired index
   **/
  SLEQP_EXPORT double* sleqp_sparse_vector_at(SleqpSparseVec* vec,
                                              int index);

  /**
   * Returns whether this vector is boxed, i.e., \f$ lb \leq x \leq ub \f$
   * for all components.
   *
   * @param[in] x    A pointer to the vector
   * @param[in] lb   A pointer to the lower bound vector
   * @param[in] ub   A pointer to the upper bound vector
   *
   SLEQP_EXPORT * @sa sleqp_sparse_vector_clip
  **/
  SLEQP_EXPORT bool sleqp_sparse_vector_is_boxed(const SleqpSparseVec* x,
                                                 const SleqpSparseVec* lb,
                                                 const SleqpSparseVec* ub);

  /**
   * Clips this vector to the specified lower and upper bounds, storing
   * the result in the pointer to a new vector, which will be boxed
   * w.r.t. `lb` and `ub`
   *
   * @param[in]  x     A pointer to the vector
   * @param[in]  lb    A pointer to the lower bound vector
   * @param[in]  ub    A pointer to the upper bound vector
   * @param[out] xclip A pointer to the clipped vector
   *
   SLEQP_EXPORT * @sa sleqp_sparse_vector_is_boxed
  **/
  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_sparse_vector_clip(const SleqpSparseVec* x,
                                         const SleqpSparseVec* lb,
                                         const SleqpSparseVec* ub,
                                         const double eps,
                                         SleqpSparseVec* xclip);

  /**
   * Prints this vector to the given file
   *
   * @param[in]  vec     A pointer to the vector
   * @param[in]  output  A pointer to an output `FILE*`
   **/
  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_sparse_vector_fprintf(const SleqpSparseVec* vec,
                                            FILE* output);

  /**
   * Returns whether the given vector is *valid*, i.e., whether
   * - `nnz` non-negative and less than or equal to `nnz_max`
   * - all indices are non-negative and less than or equal to `dim`
   * - the entries are ordered according to their indices
   * - the stored `data` is free of (IEEE) infs and NaNs
   *
   **/
  SLEQP_EXPORT bool sleqp_sparse_vector_valid(const SleqpSparseVec* vec);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_sparse_vector_free(SleqpSparseVec** vec);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_SPARSE_VEC_H */
