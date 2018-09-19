#ifndef SLEQP_FUNC_H
#define SLEQP_FUNC_H

/** @file sleqp_func.h
 * Definition of functions used for objective / constraints.
 **/

#include "sleqp_types.h"
#include "sparse/sleqp_sparse.h"

#include <assert.h>

#ifdef __cplusplus
extern "C" {
#endif

  // TODO: add SLEQPProblem* at some point

  typedef struct SleqpFunc SleqpFunc;

  /**
   * Sets the current input vector
   *
   * @param[in]  x               The input vector \f$ x \f$
   * @param[in]  num_variables   The number of variables
   * @param[out] func_grad_nnz   The number of nonzeros of the function gradient \f$ \nabla f(x) \f$
   * @param[out] cons_val_nnz    The number of nonzeros of the constraint function \f$ c(x) \f$
   * @param[out] cons_jac_nnz    The number of nonzeros of the constraint Jacobian \f$ J_c(x) \f$
   * @param[in,out] func_data    The function data
   **/
  typedef SLEQP_RETCODE (*SLEQP_FUNC_SET)(SleqpSparseVec* x,
                                          size_t num_variables,
                                          size_t* func_grad_nnz,
                                          size_t* cons_val_nnz,
                                          size_t* cons_jac_nnz,
                                          void* func_data);

  /**
   * Evaluates the function and its gradient at the current input vector
   *
   * @param[in]     num_variables   The number of variables
   * @param[out]    func_val        The function value \f$ f(x) \f$
   * @param[out]    func_grad       The function gradient \f$ \nabla f(x) \f$
   * @param[out]    cons_val        The value of the constraint function \f$ c(x) \f$
   * @param[out]    cons_jac        The constraint Jacobian \f$ J_c(x) \f$
   * @param[in,out] func_data       The function data
   **/
  typedef SLEQP_RETCODE (*SLEQP_FUNC_EVAL)(size_t num_variables,
                                           int* indices,
                                           double* func_val,
                                           SleqpSparseVec* func_grad,
                                           SleqpSparseVec* cons_val,
                                           SleqpSparseMatrix* cons_jac,
                                           void* func_data);

  /**
   * Evaluates the bilinear product of the Hessian of the Lagrangian function.
   * The Lagrangian function is given by:
   *
   * \f[
   * L(x, \lambda_0, \lambda) := \lambda_0 f(x) + \langle \lambda, c(x) \rangle
   * \f]
   *
   * The bilinear product with a direction \f$ d \f$ is then:
   * \f[
   * \frac{1}{2} \langle d, \nabla_{xx} L(x, \lambda_0, \lambda) d \rangle
   * = \frac{1}{2} \left( \lambda_0 \langle d, \nabla_{xx} f(x) d \rangle
   *   + \sum_{i=1}^{m} \lambda_i \langle d , \nabla_{xx} c_i(x) d \rangle  \right)
   * \f]
   *
   * @param[in]     num_variables     The number of variables
   * @param[in]     func_dual         The value \f$ \lambda_0 \f$
   * @param[in]     direction         The direction \f$ d \f$
   * @param[in]     cons_duals        The values \f$ \lambda \f$
   * @param[in]     bilinear_prod     The resulting bilinear product
   * @param[in,out] func_data         The function data
   *
   */
  typedef SLEQP_RETCODE (*SLEQP_HESS_EVAL_BILINEAR)(size_t num_variables,
                                                    double* func_dual,
                                                    SleqpSparseVec* direction,
                                                    SleqpSparseVec* cons_duals,
                                                    double* bilinear_prod,
                                                    void* func_data);

  /**
   * Creates a new function.
   *
   * @param[out] fstar            A pointer to the function to be created
   * @param[in]  setx             A callback to set the input vector
   * @param[in]  eval             A callback to evaluate the function and gradient
   * @param[in]  eval_bilin       A callback to evaluate the Hessian bilinear product
   * @param[in]  num_variables    The number of variables
   * @param[in]  func_data        The function data
   **/
  SLEQP_RETCODE sleqp_func_create(SleqpFunc** fstar,
                                  SLEQP_FUNC_SET setx,
                                  SLEQP_FUNC_EVAL eval,
                                  SLEQP_HESS_EVAL_BILINEAR eval_bilin,
                                  size_t num_variables,
                                  void* func_data);

  /**
   * Sets the current input vector of a function
   *
   * @param[in]  func            The function
   * @param[in]  x               The input vector \f$ x \f$
   * @param[out] func_grad_nnz   The number of nonzeros of the function gradient \f$ \nabla f(x) \f$
   * @param[out] cons_val_nnz    The number of nonzeros of the constraint function \f$ c(x) \f$
   * @param[out] cons_jac_nnz    The number of nonzeros of the constraint Jacobian \f$ J_c(x) \f$
   **/
  SLEQP_RETCODE sleqp_func_set_value(SleqpFunc* func,
                                     SleqpSparseVec* x,
                                     size_t* func_grad_nnz,
                                     size_t* cons_val_nnz,
                                     size_t* cons_jac_nnz);

  /**
   * Evaluates the function and its gradient at the current input vector
   *
   * @param[in]     func            The function
   * @param[out]    func_grad       The function gradient \f$ \nabla f(x) \f$
   * @param[out]    cons_val        The value of the constraint function \f$ c(x) \f$
   * @param[out]    cons_jac        The constraint Jacobian \f$ J_c(x) \f$
   * @param[in,out] func_data       The function data
   **/
  SLEQP_RETCODE sleqp_func_eval(SleqpFunc* func,
                                int* indices,
                                double* func_val,
                                SleqpSparseVec* func_grad,
                                SleqpSparseVec* cons_val,
                                SleqpSparseMatrix* cons_jac);

  /**
   * Evaluates the bilinear product of the Hessian of the Lagrangian of the given function.
   *
   * @param[in]     func              The function
   * @param[in]     func_dual         The value \f$ \lambda_0 \f$
   * @param[in]     direction         The direction \f$ d \f$
   * @param[in]     cons_duals        The values \f$ \lambda \f$
   * @param[in]     bilinear_prod     The resulting bilinear product
   *
   */
  SLEQP_RETCODE sleqp_hess_eval_bilinear(SleqpFunc* func,
                                         double* func_dual,
                                         SleqpSparseVec* direction,
                                         SleqpSparseVec* cons_duals,
                                         double* bilinear_prod,
                                         void* func_data);

  SLEQP_RETCODE sleqp_func_free(SleqpFunc** fstar);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_FUNC_H */
