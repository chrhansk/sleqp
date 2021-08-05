#ifndef SLEQP_FUNC_H
#define SLEQP_FUNC_H

/**
 * @file func.h
 * @brief Definition of functions used for objective / constraints.
 **/

#include "pub_func.h"

#include "timer.h"
#include "hess_struct.h"

#include <assert.h>

#ifdef __cplusplus
extern "C" {
#endif

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
  SLEQP_NODISCARD SLEQP_RETCODE sleqp_func_set_value(SleqpFunc* func,
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
  SLEQP_NODISCARD SLEQP_RETCODE sleqp_func_eval(SleqpFunc* func,
                                                const SleqpSparseVec* cons_indices,
                                                double* func_val,
                                                SleqpSparseVec* func_grad,
                                                SleqpSparseVec* cons_val,
                                                SleqpSparseMatrix* cons_jac);

  SLEQP_NODISCARD SLEQP_RETCODE sleqp_func_val(SleqpFunc* func,
                                               double* func_val);

  SLEQP_NODISCARD SLEQP_RETCODE sleqp_func_grad(SleqpFunc* func,
                                                SleqpSparseVec* func_grad);

  SLEQP_NODISCARD SLEQP_RETCODE sleqp_func_cons_val(SleqpFunc* func,
                                                    const SleqpSparseVec* cons_indices,
                                                    SleqpSparseVec* cons_val);

  SLEQP_NODISCARD SLEQP_RETCODE sleqp_func_cons_jac(SleqpFunc* func,
                                                    const SleqpSparseVec* cons_indices,
                                                    SleqpSparseMatrix* cons_jac);

  bool sleqp_func_has_psd_hessian(SleqpFunc* func);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_func_set_psd_hessian(SleqpFunc* func,
                                           bool value);

  SLEQP_FUNC_TYPE sleqp_func_get_type(SleqpFunc* func);

  SLEQP_RETCODE sleqp_func_set_type(SleqpFunc* func,
                                    SLEQP_FUNC_TYPE func_type);

  /**
   * Returns the setting timer of this function. This timer records
   * the setting of function values.
   *
   * @param[in]     func            The function
   *
   **/
  SleqpTimer* sleqp_func_get_set_timer(SleqpFunc* func);

  /**
   * Returns the evaluation timer of this function. This timer records
   * the evaluations of function values.
   *
   * @param[in]     func            The function
   *
   **/
  SleqpTimer* sleqp_func_get_val_timer(SleqpFunc* func);

  /**
   * Returns the gradient timer of this function. This timer records
   * the evaluations of function gradients.
   *
   * @param[in]     func            The function
   *
   **/
  SleqpTimer* sleqp_func_get_grad_timer(SleqpFunc* func);

  /**
   * Returns the constraint timer of this function. This timer records
   * the evaluations of constraint values.
   *
   * @param[in]     func            The function
   *
   **/
  SleqpTimer* sleqp_func_get_cons_val_timer(SleqpFunc* func);

  /**
   * Returns the Jacobian timer of this function. This timer records
   * the evaluations of constraint Jacobians.
   *
   * @param[in]     func            The function
   *
   **/
  SleqpTimer* sleqp_func_get_cons_jac_timer(SleqpFunc* func);

  /**
   * Returns the Hessian timer of this function. This timer records
   * the evaluations of Hessian products.
   *
   * @param[in]     func            The function
   *
   **/
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
  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_func_hess_prod(SleqpFunc* func,
                                     const double* func_dual,
                                     const SleqpSparseVec* direction,
                                     const SleqpSparseVec* cons_duals,
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
  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_func_hess_bilinear(SleqpFunc* func,
                                         const double* func_dual,
                                         const SleqpSparseVec* direction,
                                         const SleqpSparseVec* cons_duals,
                                         double* bilinear_prod);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_FUNC_H */
