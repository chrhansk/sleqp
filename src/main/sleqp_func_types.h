#ifndef SLEQP_FUNC_TYPES_H
#define SLEQP_FUNC_TYPES_H

#include "sleqp_types.h"

#include "sparse/sleqp_sparse_matrix.h"
#include "sparse/sleqp_sparse_vec.h"

/**
 * @file sleqp_func_types.h
 * @brief Definition of function types used for objective / constraints.
 **/

#ifdef __cplusplus
extern "C" {
#endif

  typedef enum {
    SLEQP_VALUE_REASON_NONE,
    SLEQP_VALUE_REASON_INIT,
    SLEQP_VALUE_REASON_CHECKING_DERIV,
    SLEQP_VALUE_REASON_ACCEPTED_ITERATE,
    SLEQP_VALUE_REASON_TRYING_ITERATE,
    SLEQP_VALUE_REASON_TRYING_SOC_ITERATE,
    SLEQP_VALUE_REASON_REJECTED_ITERATE,
  } SLEQP_VALUE_REASON;

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
                                          SLEQP_VALUE_REASON reason,
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
   **/
  typedef SLEQP_RETCODE (*SLEQP_HESS_PRODUCT)(int num_variables,
                                              double* func_dual,
                                              SleqpSparseVec* direction,
                                              SleqpSparseVec* cons_duals,
                                              SleqpSparseVec* product,
                                              void* func_data);

  /**
   * Cleans up any allocated memory stored in the function data.
   *
   * @param[in,out] func_data  The function data
   *
   **/
  typedef SLEQP_RETCODE (*SLEQP_FUNC_FREE)(void* func_data);

  typedef struct {
    SLEQP_FUNC_SET set_value;
    SLEQP_FUNC_EVAL func_eval;
    SLEQP_HESS_PRODUCT hess_prod;
    SLEQP_FUNC_FREE func_free;
  } SleqpFuncCallbacks;

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_FUNC_TYPES_H */
