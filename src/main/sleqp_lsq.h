#ifndef SLEQP_LSQ_H
#define SLEQP_LSQ_H

/**
 * @file sleqp_lsq.h
 * @brief Defintion of least squares functions.
 **/

#include "sleqp_export.h"
#include "sleqp_iterate.h"
#include "sleqp_problem.h"
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
   * @param[in]     func              The function
   * @param[out]    residual          The resulting residual
   * @param[in,out] func_data         The function data
   *
   */
  typedef SLEQP_RETCODE (*SLEQP_LSQ_EVAL)(SleqpFunc* func,
                                          SleqpSparseVec* residual,
                                          void* func_data);

  /**
   * Evaluates the product of the Jacobian of the constraints part and a direction.
   * \f[
   *    \nabla_{x} c(x) d
   * \f]
   * This is useful when the constraints part is used to represent a least-squares residual.
   *
   * @param[in]     func              The function
   * @param[in]     forward_direction The direction \f$ d \f$
   * @param[out]    product           The resulting product
   * @param[in,out] func_data         The function data
   *
   */
  typedef SLEQP_RETCODE (*SLEQP_LSQ_JAC_FORWARD)(SleqpFunc* func,
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
   * @param[in]     func              The function
   * @param[in]     adjoint_direction The direction \f$ d \f$
   * @param[out]    product           The resulting product
   * @param[in,out] func_data         The function data
   *
   */
  typedef SLEQP_RETCODE (*SLEQP_LSQ_JAC_ADJOINT)(SleqpFunc* func,
                                                 const SleqpSparseVec* adjoint_direction,
                                                 SleqpSparseVec* product,
                                                 void* func_data);

  typedef struct {
    SLEQP_FUNC_SET set_value;
    SLEQP_LSQ_EVAL lsq_eval;
    SLEQP_LSQ_JAC_FORWARD lsq_jac_forward;
    SLEQP_LSQ_JAC_ADJOINT lsq_jac_adjoint;
    SLEQP_FUNC_VAL additional_func_val;
    SLEQP_FUNC_GRAD additional_func_grad;
    SLEQP_FUNC_CONS_VAL additional_cons_val;
    SLEQP_FUNC_CONS_JAC additional_cons_jac;
    SLEQP_HESS_PROD additional_hess_prod;
    SLEQP_FUNC_FREE func_free;
  } SleqpLSQCallbacks;

  SLEQP_EXPORT SLEQP_RETCODE sleqp_lsq_func_create(SleqpFunc** fstar,
                                                   SleqpLSQCallbacks* callbacks,
                                                   int num_variables,
                                                   int num_constraints,
                                                   int num_residuals,
                                                   double levenberg_marquardt,
                                                   SleqpParams* params,
                                                   void* func_data);

  /**
   * Sets the callbacks of this LSQ function to the specified ones
   **/
  SLEQP_EXPORT SLEQP_RETCODE sleqp_lsq_func_set_callbacks(SleqpFunc* func,
                                                          SleqpLSQCallbacks* callbacks);


#ifdef __cplusplus
}
#endif

#endif /* SLEQP_LSQ_H */
