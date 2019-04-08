#ifndef SLEQP_FUNC_H
#define SLEQP_FUNC_H

/**
 * @file sleqp_func.h
 * @brief Definition of functions used for objective / constraints.
 **/

#include "sleqp_types.h"

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
                                          int num_variables,
                                          int* func_grad_nnz,
                                          int* cons_val_nnz,
                                          int* cons_jac_nnz,
                                          void* func_data);

  /**
   * Evaluates the function and its gradient at the current input vector
   *
   * @param[in]     num_variables   The number of variables
   * @param[in]     cons_indices    The indices of the constraint function
   *                                to be evaluated
   * @param[out]    func_val        The function value \f$ f(x) \f$
   * @param[out]    func_grad       The function gradient \f$ \nabla f(x) \f$
   * @param[out]    cons_val        The value of the constraint function \f$ c(x) \f$
   * @param[out]    cons_jac        The constraint Jacobian \f$ J_c(x) \f$
   * @param[in,out] func_data       The function data
   **/
  typedef SLEQP_RETCODE (*SLEQP_FUNC_EVAL)(int num_variables,
                                           SleqpSparseVec* cons_indices,
                                           double* func_val,
                                           SleqpSparseVec* func_grad,
                                           SleqpSparseVec* cons_val,
                                           SleqpSparseMatrix* cons_jac,
                                           void* func_data);

  /**
   * Evaluates the product of the Hessian of the Lagrangian function.
   * The Lagrangian function is given by:
   *
   * \f[
   * L(x, \lambda_0, \lambda) := \lambda_0 f(x) + \langle \lambda, c(x) \rangle
   * \f]
   *
   * The product with a direction \f$ d \f$ is then:
   * \f[
   * \nabla_{xx} L(x, \lambda_0, \lambda) d
   * = \left( \lambda_0 \nabla_{xx} f(x) d
   *   + \sum_{i=1}^{m} \lambda_i  \nabla_{xx} c_i(x) d \right)
   * \f]
   *
   * @param[in]     num_variables     The number of variables
   * @param[in]     func_dual         The value \f$ \lambda_0 \f$
   * @param[in]     direction         The direction \f$ d \f$
   * @param[in]     cons_duals        The values \f$ \lambda \f$
   * @param[out]    product           The resulting product
   * @param[in,out] func_data         The function data
   *
   */
  typedef SLEQP_RETCODE (*SLEQP_HESS_PRODUCT)(int num_variables,
                                              double* func_dual,
                                              SleqpSparseVec* direction,
                                              SleqpSparseVec* cons_duals,
                                              SleqpSparseVec* product,
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
                                  SLEQP_HESS_PRODUCT eval_hess_prod,
                                  int num_variables,
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

  /**
   * Returns the number of variables \f$ n \f$.
   **/
  int sleqp_func_get_num_variables(SleqpFunc* func);

  /**
   * Returns the number of evaluations of the function \f$ f \f$.
   **/
  int sleqp_func_get_num_func_evals(SleqpFunc* func);

  /**
   * Returns the number of evaluations of the function \f$ c \f$.
   **/
  int sleqp_func_get_num_cons_evals(SleqpFunc* func);

  /**
   * Returns the number of evaluations of the gradient of \f$ f \f$.
   **/
  int sleqp_func_get_num_grad_evals(SleqpFunc* func);

  /**
   * Returns the number of evaluations of the Jacobian of \f$ c \f$.
   **/
  int sleqp_func_get_num_jac_evals(SleqpFunc* func);

  /**
   * Returns the number of evaluations of Hessian products of the Lagrangian \f$ L \f$.
   *
   **/
  int sleqp_func_get_num_hess_evals(SleqpFunc* func);

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

  /**
   * Frees a previously created function.
   **/
  SLEQP_RETCODE sleqp_func_free(SleqpFunc** fstar);

  /**
   * @}
   **/

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_FUNC_H */
