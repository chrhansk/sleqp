#ifndef SLEQP_HESS_STRUCT_H
#define SLEQP_HESS_STRUCT_H

#include "sleqp_export.h"
#include "sleqp_types.h"

/**
 * @file sleqp_solver.h
 * @brief Definition of problem block structure.
 **/

#ifdef __cplusplus
extern "C" {
#endif
  /**
   *
   * Hessian structure of the Lagrangian
   *
   * \f[ L(x, \lambda, \mu) = f(x) + \langle \lambda, c \rangle + \langle 1, \mu \rangle \f].
   *
   * The Hessian \f$ H_L \f$ is assumed to consist of a number of \f$k\f$ blocks,
   * given in terms of indices \f$ 1 = j_1 < j_2 < \ldots < j_{k+1} \leq n \f$.
   * All non-zero entries \f$ (i, j) \f$ must satisfy that
   * \f$ j_l \leq i, j < j_{l + 1} \f$ for some \f$ l = 1, \ldots, k \f$.
   *
   * The default is to assume one block of size \f$ n \f$, implying
   * no particular structure of the Hessian.
   *
   * If \f$ j_{k+1} < n \f$, then the range from \f$ j_{k+1} \f$ to \f$ n \f$
   * must only contain zero entries, i.e., combinations of linear constraints
   * and variables.
   *
   **/
  typedef struct SleqpHessianStruct SleqpHessianStruct;

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_hessian_struct_create(SleqpHessianStruct** star,
                                            int dimension,
                                            bool empty);

  /**
   * Returns the number \f$ k \f$ of blocks.
   *
   * @param[in]  hessian_struct  The Hessian structure
   * @returns                    The number \f$ k \f$ of blocks
   **/
  SLEQP_EXPORT int sleqp_hessian_struct_get_num_blocks(const SleqpHessianStruct* hessian_struct);

  /**
   * Returns the \f$ l \f$-th block of the Hessian
   *
   * @param[in]  hessian_struct  The Hessian structure
   * @param[in]  block           The index \f$ l \f$ of the block
   * @param[out] begin           The 0-based index \f$ j_l \f$
   * @param[out] end             The 0-based index \f$ j_{l+1} \f$
   **/
  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_hessian_struct_get_block_range(const SleqpHessianStruct* hessian_struct,
                                                     int block,
                                                     int* begin,
                                                     int* end);

  /**
   * Pushes a new block into the Hessian
   *
   * @param[in]  hessian_struct  The Hessian structure
   * @param[out] end             The 0-based index \f$ j_{l+1} \f$
   **/
  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_hessian_struct_push_block(SleqpHessianStruct* hessian_struct,
                                                int end);

  /**
   * Clears the Hessian structure, i.e., sets \f$ k = 0 \f$
   *
   * @param[in]  hessian_struct  The Hessian structure
   **/
  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_hessian_struct_clear(SleqpHessianStruct* hessian_struct);

  /**
   * Returns the linear range
   *
   * @param[in]   hessian_struct  The Hessian structure
   * @param[out]  begin           The value \f$ j_{k + 1} \f$
   * @param[out]  end             The value \f$ n \f$
   *
   **/
  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_hessian_struct_get_linear_range(const SleqpHessianStruct* hessian_struct,
                                                      int* begin,
                                                      int* end);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_hessian_struct_copy(const SleqpHessianStruct* source,
                                          SleqpHessianStruct* target);

  /**
   * Prints the Hessian structure
   **/
  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_hessian_struct_fprintf(SleqpHessianStruct* hessian_struct,
                                             FILE* output);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_hessian_struct_capture(SleqpHessianStruct* hessian_struct);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_hessian_struct_release(SleqpHessianStruct** star);


#ifdef __cplusplus
}
#endif

#endif /* SLEQP_HESS_STRUCT_H */
