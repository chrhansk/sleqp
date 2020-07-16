#ifndef SLEQP_FUNC_H
#define SLEQP_FUNC_H

/**
 * @file sleqp_func.h
 * @brief Definition of functions used for objective / constraints.
 **/

#include "sleqp_types.h"
#include "sleqp_func_types.h"

#include "sleqp_timer.h"
#include "sleqp_hess_struct.h"

#include "sparse/sleqp_sparse_matrix.h"
#include "sparse/sleqp_sparse_vec.h"

#include <assert.h>

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpFunc SleqpFunc;

  /**
   * @defgroup function Function definition
   * @{
   **/

  /**
   * Creates a new function.
   *
   * @param[out] fstar            A pointer to the function to be created
   * @param[in]  setx             A callback to set the input vector
   * @param[in]  eval             A callback to evaluate the function and gradient
   * @param[in]  num_variables    The number of variables
   * @param[in]  func_data        The function data
   **/
  SLEQP_RETCODE sleqp_func_create(SleqpFunc** fstar,
                                  SleqpFuncCallbacks* callbacks,
                                  int num_variables,
                                  void* func_data);

  /**
   * Sets the current input vector of a function
   *
   * @param[in]  func            The function
   * @param[in]  x               The input vector \f$ x \f$
   * @param[in]  reason          The reason for setting \f$ x \f$
   * @param[out] func_grad_nnz   The number of nonzeros of the function gradient \f$ \nabla f(x) \f$
   * @param[out] cons_val_nnz    The number of nonzeros of the constraint function \f$ c(x) \f$
   * @param[out] cons_jac_nnz    The number of nonzeros of the constraint Jacobian \f$ J_c(x) \f$
   **/
  SLEQP_RETCODE sleqp_func_set_value(SleqpFunc* func,
                                     SleqpSparseVec* x,
                                     SLEQP_VALUE_REASON reason,
                                     int* func_grad_nnz,
                                     int* cons_val_nnz,
                                     int* cons_jac_nnz);

  /**
   * Evaluates the function and its gradient at the current input vector
   *
   * @param[in]     func            The function
   * @param[in]     cons_indices    The indices of the constraint function
   *                                to be evaluated
   * @param[out]    func_grad       The function gradient \f$ \nabla f(x) \f$
   * @param[out]    cons_val        The value of the constraint function \f$ c(x) \f$
   * @param[out]    cons_jac        The constraint Jacobian \f$ J_c(x) \f$
   * @param[in,out] func_data       The function data
   **/
  SLEQP_RETCODE sleqp_func_eval(SleqpFunc* func,
                                SleqpSparseVec* cons_indices,
                                double* func_val,
                                SleqpSparseVec* func_grad,
                                SleqpSparseVec* cons_val,
                                SleqpSparseMatrix* cons_jac);

  SleqpHessianStruct* sleqp_func_get_hess_struct(SleqpFunc* func);

  /**
   * Returns the number of variables \f$ n \f$.
   **/
  int sleqp_func_get_num_variables(SleqpFunc* func);

  /**
   * Returns the number of evaluations of Hessian products of the Lagrangian \f$ L \f$.
   *
   **/
  int sleqp_func_get_num_hess_evals(SleqpFunc* func);

  SleqpTimer* sleqp_func_get_eval_timer(SleqpFunc* func);

  SleqpTimer* sleqp_func_get_hess_timer(SleqpFunc* func);

  /**
   * Evaluates the product of the Hessian of the Lagrangian of the given function.
   *
   * @param[in]     func              The function
   * @param[in]     func_dual         The value \f$ \lambda_0 \f$
   * @param[in]     direction         The direction \f$ d \f$
   * @param[in]     cons_duals        The values \f$ \lambda \f$
   * @param[out]    product           The resulting product
   *
   */
  SLEQP_RETCODE sleqp_func_hess_prod(SleqpFunc* func,
                                     double* func_dual,
                                     SleqpSparseVec* direction,
                                     SleqpSparseVec* cons_duals,
                                     SleqpSparseVec* product);

  /**
   * Evaluates the bilinear product of the Hessian of the Lagrangian of the given function.
   *
   * @param[in]     func              The function
   * @param[in]     func_dual         The value \f$ \lambda_0 \f$
   * @param[in]     direction         The direction \f$ d \f$
   * @param[in]     cons_duals        The values \f$ \lambda \f$
   * @param[out]    bilinear_prod     The resulting bilinear product
   *
   */
  SLEQP_RETCODE sleqp_func_hess_bilinear(SleqpFunc* func,
                                         double* func_dual,
                                         SleqpSparseVec* direction,
                                         SleqpSparseVec* cons_duals,
                                         double* bilinear_prod);

  /**
   * Returns the function data associated with the given function.
   **/
  void* sleqp_func_get_data(SleqpFunc* func);

  SLEQP_RETCODE sleqp_func_capture(SleqpFunc* func);

  SLEQP_RETCODE sleqp_func_release(SleqpFunc** fstar);

  /**
   * @}
   **/

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_FUNC_H */
