#ifndef SLEQP_LSQ_H
#define SLEQP_LSQ_H

/**
 * @file sleqp_lsq.h
 * @brief Defintion of least squares functions.
 **/

#include "sleqp_problem.h"

#include "sleqp_iterate.h"
#include "sparse/sleqp_sparse_vec.h"

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Evaluates the residual.
   * \f[
   *    \nabla_{x} c(x) d
   * \f]
   * This is useful when the constraints part is used to represent a least-squares residual.
   *
   * @param[in]     num_variables     The number of variables
   * @param[out]    residual          The resulting residual
   * @param[in,out] func_data         The function data
   *
   */
  typedef SLEQP_RETCODE (*SLEQP_LSQ_EVAL)(int num_variables,
                                          SleqpSparseVec* residual,
                                          void* func_data);

  /**
   * Evaluates the product of the Jacobian of the constraints part and a direction.
   * \f[
   *    \nabla_{x} c(x) d
   * \f]
   * This is useful when the constraints part is used to represent a least-squares residual.
   *
   * @param[in]     num_variables     The number of variables
   * @param[in]     forward_direction The direction \f$ d \f$
   * @param[out]    product           The resulting product
   * @param[in,out] func_data         The function data
   *
   */
  typedef SLEQP_RETCODE (*SLEQP_LSQ_JAC_FORWARD)(int num_variables,
                                                 const SleqpSparseVec* forward_direction,
                                                 SleqpSparseVec* product,
                                                 void* func_data);

  /**
   * Evaluates the product of the transposed Jacobian of the constraints part and a direction.
   * \f[
   *    d^T \nabla_{x} c(x)
   * \f]
   * This is useful when the constraints part is used to represent a least-squares residual.
   *
   * @param[in]     num_variables     The number of variables
   * @param[in]     adjoint_direction The direction \f$ d \f$
   * @param[out]    product           The resulting product
   * @param[in,out] func_data         The function data
   *
   */
  typedef SLEQP_RETCODE (*SLEQP_LSQ_JAC_ADJOINT)(int num_variables,
                                                 const SleqpSparseVec* adjoint_direction,
                                                 SleqpSparseVec* product,
                                                 void* func_data);

  typedef struct {
    SLEQP_FUNC_SET set_value;
    SLEQP_LSQ_EVAL lsq_eval;
    SLEQP_LSQ_JAC_FORWARD lsq_jac_forward;
    SLEQP_LSQ_JAC_ADJOINT lsq_jac_adjoint;
    SLEQP_FUNC_EVAL eval_additional;
    SLEQP_HESS_PRODUCT hess_prod_additional;
    SLEQP_FUNC_FREE func_free;
  } SleqpLSQCallbacks;

  SLEQP_RETCODE sleqp_lsq_func_create(SleqpFunc** fstar,
                                      SleqpLSQCallbacks* callbacks,
                                      int num_variables,
                                      int num_residuals,
                                      double levenberg_marquardt,
                                      SleqpParams* params,
                                      void* func_data);

  /**
   * Sets the callbacks of this LSQ function to the specified ones
   **/
  SLEQP_RETCODE sleqp_lsq_func_set_callbacks(SleqpFunc* func,
                                             SleqpLSQCallbacks* callbacks);


#ifdef __cplusplus
}
#endif

#endif /* SLEQP_LSQ_H */
