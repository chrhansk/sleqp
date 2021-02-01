#ifndef SLEQP_FUNC_H
#define SLEQP_FUNC_H

/**
 * @file sleqp_func.h
 * @brief Definition of functions used for objective / constraints.
 **/

#include "sleqp_export.h"
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
   * Creates a new function
   *
   * @param[out] fstar            A pointer to the function to be created
   * @param[in]  setx             A callback to set the input vector
   * @param[in]  eval             A callback to evaluate the function and gradient
   * @param[in]  num_variables    The number of variables
   * @param[in]  num_constraints  The number of constraints
   * @param[in]  func_data        The function data
   **/
  SLEQP_EXPORT SLEQP_RETCODE sleqp_func_create(SleqpFunc** fstar,
                                               SleqpFuncCallbacks* callbacks,
                                               int num_variables,
                                               int num_constraints,
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
  SLEQP_EXPORT SLEQP_RETCODE sleqp_func_set_value(SleqpFunc* func,
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
  SLEQP_EXPORT SLEQP_RETCODE sleqp_func_eval(SleqpFunc* func,
                                             const SleqpSparseVec* cons_indices,
                                             double* func_val,
                                             SleqpSparseVec* func_grad,
                                             SleqpSparseVec* cons_val,
                                             SleqpSparseMatrix* cons_jac);

  SLEQP_EXPORT SLEQP_RETCODE sleqp_func_val(SleqpFunc* func,
                                             double* func_val);

  SLEQP_EXPORT SLEQP_RETCODE sleqp_func_grad(SleqpFunc* func,
                                             SleqpSparseVec* func_grad);

  SLEQP_EXPORT SLEQP_RETCODE sleqp_func_cons_val(SleqpFunc* func,
                                                 const SleqpSparseVec* cons_indices,
                                                 SleqpSparseVec* cons_val);

  SLEQP_EXPORT SLEQP_RETCODE sleqp_func_cons_jac(SleqpFunc* func,
                                                 const SleqpSparseVec* cons_indices,
                                                 SleqpSparseMatrix* cons_jac);

  /**
   * Sets the callbacks of this function to the specified ones
   **/
  SLEQP_EXPORT SLEQP_RETCODE sleqp_func_set_callbacks(SleqpFunc* func,
                                                      SleqpFuncCallbacks* callbacks);

  /**
   * Returns the Hessian structure of this function
   *
   * @param[in]     func            The function
   *
   **/
  SLEQP_EXPORT SleqpHessianStruct* sleqp_func_get_hess_struct(SleqpFunc* func);

  /**
   * Returns the number of variables \f$ n \f$.
   **/
  SLEQP_EXPORT int sleqp_func_get_num_variables(SleqpFunc* func);

  /**
   * Returns the number of constraints \f$ m \f$.
   **/
  SLEQP_EXPORT int sleqp_func_get_num_constraints(SleqpFunc* func);

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
  SLEQP_EXPORT SleqpTimer* sleqp_func_get_hess_timer(SleqpFunc* func);

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
  SLEQP_EXPORT SLEQP_RETCODE sleqp_func_hess_prod(SleqpFunc* func,
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
  SLEQP_EXPORT SLEQP_RETCODE sleqp_func_hess_bilinear(SleqpFunc* func,
                                                      const double* func_dual,
                                                      const SleqpSparseVec* direction,
                                                      const SleqpSparseVec* cons_duals,
                                                      double* bilinear_prod);

  /**
   * Returns the function data associated with the given function.
   **/
  SLEQP_EXPORT void* sleqp_func_get_data(SleqpFunc* func);

  SLEQP_EXPORT SLEQP_RETCODE sleqp_func_capture(SleqpFunc* func);

  SLEQP_EXPORT SLEQP_RETCODE sleqp_func_release(SleqpFunc** fstar);

  /**
   * @}
   **/

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_FUNC_H */
