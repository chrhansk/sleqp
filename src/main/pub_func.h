#ifndef PUB_FUNC_H
#define PUB_FUNC_H

#include "pub_hess_struct.h"
#include "sparse/pub_sparse_matrix.h"
#include "sparse/pub_sparse_vec.h"

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * @defgroup function Function definition
   * @{
   **/

  typedef struct SleqpFunc SleqpFunc;

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
   * @param[in]     func            The function
   * @param[in]     value           The value
   * @param[out]    func_grad_nnz   The number of nonzeros of the function gradient \f$ \nabla f(x) \f$
   * @param[out]    cons_val_nnz    The number of nonzeros of the constraint function \f$ c(x) \f$
   * @param[out]    cons_jac_nnz    The number of nonzeros of the constraint Jacobian \f$ J_c(x) \f$
   * @param[in,out] func_data       The function data
   **/
  typedef SLEQP_RETCODE (*SLEQP_FUNC_SET)(SleqpFunc* func,
                                          SleqpSparseVec* value,
                                          SLEQP_VALUE_REASON reason,
                                          int* func_grad_nnz,
                                          int* cons_val_nnz,
                                          int* cons_jac_nnz,
                                          void* func_data);

  /**
   * Evaluates the function and its gradient at the current input vector
   *
   * @param[in]     func            The function
   * @param[in]     cons_indices    The indices of the constraint function
   *                                to be evaluated
   * @param[out]    func_val        The function value \f$ f(x) \f$
   * @param[out]    func_grad       The function gradient \f$ \nabla f(x) \f$
   * @param[out]    cons_val        The value of the constraint function \f$ c(x) \f$
   * @param[out]    cons_jac        The constraint Jacobian \f$ J_c(x) \f$
   * @param[in,out] func_data       The function data
   **/

  /**
   * Evaluates the function at the current input vector
   *
   * @param[in]     func            The function
   * @param[out]    func_val        The function value \f$ f(x) \f$
   * @param[in,out] func_data       The function data
   **/
  typedef SLEQP_RETCODE (*SLEQP_FUNC_VAL)(SleqpFunc* func,
                                          double* func_val,
                                          void* func_data);

  /**
   * Evaluates the function gradient at the current input vector
   *
   * @param[in]     func            The function
   * @param[out]    func_grad       The function gradient \f$ \nabla f(x) \f$
   * @param[in,out] func_data       The function data
   **/
  typedef SLEQP_RETCODE (*SLEQP_FUNC_GRAD)(SleqpFunc* func,
                                           SleqpSparseVec* func_grad,
                                           void* func_data);

  /**
   * Evaluates the constraints at the current input vector
   *
   * @param[in]     func            The function
   * @param[in]     cons_indices    The indices of the constraint function
   *                                to be evaluated
   * @param[out]    cons_val        The value of the constraint function \f$ c(x) \f$
   * @param[in,out] func_data       The function data
   **/
  typedef SLEQP_RETCODE (*SLEQP_FUNC_CONS_VAL)(SleqpFunc* func,
                                               const SleqpSparseVec* cons_indices,
                                               SleqpSparseVec* cons_val,
                                               void* func_data);

  /**
   * Evaluates the constraing Jacobian at the current input vector
   *
   * @param[in]     func            The function
   * @param[in]     cons_indices    The indices of the constraint function
   *                                to be evaluated
   * @param[out]    cons_jac        The constraint Jacobian \f$ J_c(x) \f$
   * @param[in,out] func_data       The function data
   **/
  typedef SLEQP_RETCODE (*SLEQP_FUNC_CONS_JAC)(SleqpFunc* func,
                                               const SleqpSparseVec* cons_indices,
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
   * @param[in]     func              The function
   * @param[in]     func_dual         The value \f$ \lambda_0 \f$
   * @param[in]     direction         The direction \f$ d \f$
   * @param[in]     cons_duals        The values \f$ \lambda \f$
   * @param[out]    product           The resulting product
   * @param[in,out] func_data         The function data
   *
   **/
  typedef SLEQP_RETCODE (*SLEQP_HESS_PROD)(SleqpFunc* func,
                                           const double* func_dual,
                                           const SleqpSparseVec* direction,
                                           const SleqpSparseVec* cons_duals,
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
    SLEQP_FUNC_VAL func_val;
    SLEQP_FUNC_GRAD func_grad;
    SLEQP_FUNC_CONS_VAL cons_val;
    SLEQP_FUNC_CONS_JAC cons_jac;
    SLEQP_HESS_PROD hess_prod;
    SLEQP_FUNC_FREE func_free;
  } SleqpFuncCallbacks;

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
  SLEQP_EXPORT SLEQP_NODISCARD SLEQP_RETCODE sleqp_func_create(SleqpFunc** fstar,
                                                               SleqpFuncCallbacks* callbacks,
                                                               int num_variables,
                                                               int num_constraints,
                                                               void* func_data);

  /**
   * Returns the number of variables \f$ n \f$.
   **/
  SLEQP_EXPORT int sleqp_func_get_num_variables(SleqpFunc* func);

  /**
   * Returns the number of constraints \f$ m \f$.
   **/
  SLEQP_EXPORT int sleqp_func_get_num_constraints(SleqpFunc* func);

  /**
   * Sets the callbacks of this function to the specified ones
   **/
  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_func_set_callbacks(SleqpFunc* func,
                                         SleqpFuncCallbacks* callbacks);

  /**
   * Returns the Hessian structure of this function
   *
   * @param[in]     func            The function
   *
   **/
  SLEQP_EXPORT
  SleqpHessianStruct* sleqp_func_get_hess_struct(SleqpFunc* func);

  /**
   * Returns the function data associated with the given function.
   **/
  SLEQP_EXPORT void* sleqp_func_get_data(SleqpFunc* func);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_func_capture(SleqpFunc* func);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_func_release(SleqpFunc** fstar);

  /**
   * @}
   **/


#ifdef __cplusplus
}
#endif

#endif /* PUB_FUNC_H */