#ifndef SLEQP_PUB_LSQ_H
#define SLEQP_PUB_LSQ_H

/**
 * @file lsq.h
 * @brief Defintion of least squares functions.
 **/

#include "sleqp/export.h"
#include "sleqp/pub_iterate.h"
#include "sleqp/pub_problem.h"
#include "sleqp/sparse/pub_vec.h"

/**
 * A least-squares function consists of a residual \f$ r : \mathbb{R}^n \to
 *\mathbb{R}^k \f$ defining the objective as \f$ 1/2 \|r(x)\|^2 \f$ subjecting
 *to the usual constraints.
 **/

/**
 * Queries the number of nonzeros of the function at the
 * current primal point.
 *
 * @param[in]     func            The function
 * @param[out]    residual_nnz    The number of nonzeros of the residuals
 *\f$r(x)\f$
 * @param[out]    jac_fwd_nnz     The number of nonzeros of
 *\f$J_r(x)d\f$
 * @param[out]    jac_adj_nnz     The number of nonzeros of
 *\f$d^{T}J_r(x)\f$
 * @param[out]    cons_val_nnz    The number of nonzeros of the constraint
 *function \f$ c(x) \f$
 * @param[out]    cons_jac_nnz    The number of nonzeros of the constraint
 *Jacobian \f$ J_c(x) \f$
 * @param[out]    cons_jac_nnz    The number of nonzeros of Hessian products
 **/
typedef SLEQP_RETCODE (*SLEQP_LSQ_NONZEROS)(SleqpFunc* func,
                                            int* residual_nnz,
                                            int* jac_fwd_nnz,
                                            int* jac_adj_nnz,
                                            int* cons_val_nnz,
                                            int* cons_jac_nnz,
                                            void* func_data);

/**
 * Evaluates the residual.
 * \f[
 *    r(x)
 * \f]
 * This is useful when the constraints part is used to represent a least-squares
 * residual.
 *
 * @param[in]     func              The function
 * @param[out]    residual          The resulting residual
 * @param[in,out] func_data         The function data
 *
 */
typedef SLEQP_RETCODE (*SLEQP_LSQ_RESIDUALS)(SleqpFunc* func,
                                             SleqpVec* residual,
                                             void* func_data);

/**
 * Evaluates the product of the Jacobian of the residual with a direction \f$ d
 * \in \mathbb{R}^n \f$. \f[ J_r(x) d \f]
 *
 * @param[in]     func              The function
 * @param[in]     forward_direction The direction \f$ d \f$
 * @param[out]    product           The resulting product
 * @param[in,out] func_data         The function data
 *
 */
typedef SLEQP_RETCODE (*SLEQP_LSQ_JAC_FORWARD)(
  SleqpFunc* func,
  const SleqpVec* forward_direction,
  SleqpVec* product,
  void* func_data);

/**
 * Evaluates the product of the Jacobian of the residual with a direction \f$ d
 * \in \mathbb{R}^k \f$. \f[ d^T J_r(x) \f]
 *
 * @param[in]     func              The function
 * @param[in]     adjoint_direction The direction \f$ d \f$
 * @param[out]    product           The resulting product
 * @param[in,out] func_data         The function data
 *
 */
typedef SLEQP_RETCODE (*SLEQP_LSQ_JAC_ADJOINT)(
  SleqpFunc* func,
  const SleqpVec* adjoint_direction,
  SleqpVec* product,
  void* func_data);

typedef struct
{
  SLEQP_FUNC_SET set_value;
  SLEQP_LSQ_NONZEROS lsq_nonzeros;
  SLEQP_LSQ_RESIDUALS lsq_residuals;
  SLEQP_LSQ_JAC_FORWARD lsq_jac_forward;
  SLEQP_LSQ_JAC_ADJOINT lsq_jac_adjoint;
  SLEQP_FUNC_CONS_VAL cons_val;
  SLEQP_FUNC_CONS_JAC cons_jac;
  SLEQP_FUNC_FREE func_free;
} SleqpLSQCallbacks;

SLEQP_EXPORT SLEQP_NODISCARD SLEQP_RETCODE
sleqp_lsq_func_create(SleqpFunc** fstar,
                      SleqpLSQCallbacks* callbacks,
                      int num_variables,
                      int num_constraints,
                      int num_residuals,
                      double lm_factor,
                      SleqpParams* params,
                      void* func_data);

/**
 * Sets the callbacks of this LSQ function to the specified ones
 **/
SLEQP_EXPORT SLEQP_NODISCARD SLEQP_RETCODE
sleqp_lsq_func_set_callbacks(SleqpFunc* func, SleqpLSQCallbacks* callbacks);

SLEQP_EXPORT SLEQP_NODISCARD SLEQP_RETCODE
sleqp_lsq_func_set_lm_factor(SleqpFunc* func, double lm_factor);

#endif /* SLEQP_PUB_LSQ_H */
