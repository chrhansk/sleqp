#ifndef SLEQP_FUNC_H
#define SLEQP_FUNC_H

/**
 * @file func.h
 * @brief Definition of functions used for objective / constraints.
 **/

#include "pub_func.h"

#include "hess_struct.h"
#include "timer.h"

#include <assert.h>

typedef enum
{
  SLEQP_FUNC_HESS_INEXACT  = (1 << 0),
  SLEQP_FUNC_HESS_PSD      = (1 << 1),
  SLEQP_FUNC_INTERNAL      = (1 << 2),
  SLEQP_FUNC_HESS_INTERNAL = (1 << 3)
} SLEQP_FUNC_FLAGS;

#define SLEQP_FUNC_CALL(x, noraise, message)                                   \
  do                                                                           \
  {                                                                            \
    if ((noraise))                                                             \
    {                                                                          \
      SLEQP_CALL(x);                                                           \
    }                                                                          \
    else                                                                       \
    {                                                                          \
      SLEQP_RETCODE retcode = (x);                                             \
      if (retcode != SLEQP_OKAY)                                               \
      {                                                                        \
        sleqp_raise(SLEQP_FUNC_EVAL_ERROR, message);                           \
      }                                                                        \
    }                                                                          \
  } while (false)

#define SLEQP_FUNC_ERROR_SET_VALUE "Error setting function value"
#define SLEQP_FUNC_ERROR_OBJ_VAL "Error evaluating objective"
#define SLEQP_FUNC_ERROR_OBJ_GRAD "Error evaluating objective gradient"
#define SLEQP_FUNC_ERROR_CONS_VAL "Error evaluating constraints"
#define SLEQP_FUNC_ERROR_CONS_JAC "Error evaluating constraint Jacobian"
#define SLEQP_FUNC_ERROR_HESS_PROD "Error evaluating Hessian product"

/**
 * Sets the current input vector of a function
 *
 * @param[in]  func            The function
 * @param[in]  x               The input vector \f$ x \f$
 * @param[in]  reason          The reason for setting \f$ x \f$
 * @param[out] reject          Whether to manually reject the step
 * @param[out] obj_grad_nnz    The number of nonzeros of the function gradient
 *\f$ \nabla f(x) \f$
 * @param[out] cons_val_nnz    The number of nonzeros of the constraint function
 *\f$ c(x) \f$
 * @param[out] cons_jac_nnz    The number of nonzeros of the constraint Jacobian
 *\f$ J_c(x) \f$
 **/
SLEQP_NODISCARD SLEQP_RETCODE
sleqp_func_set_value(SleqpFunc* func,
                     SleqpVec* x,
                     SLEQP_VALUE_REASON reason,
                     bool* reject,
                     int* obj_grad_nnz,
                     int* cons_val_nnz,
                     int* cons_jac_nnz);

/**
 * Evaluates the function and its gradient at the current input vector
 *
 * @param[in]     func            The function
 * @param[out]    obj_grad        The objective gradient \f$ \nabla f(x) \f$
 * @param[out]    cons_val        The value of the constraint function \f$ c(x)
 *\f$
 * @param[out]    cons_jac        The constraint Jacobian \f$ J_c(x) \f$
 * @param[in,out] func_data       The function data
 **/
SLEQP_NODISCARD SLEQP_RETCODE
sleqp_func_eval(SleqpFunc* func,
                double* obj,
                SleqpVec* obj_grad,
                SleqpVec* cons_val,
                SleqpSparseMatrix* cons_jac);

SLEQP_NODISCARD SLEQP_RETCODE
sleqp_func_obj_val(SleqpFunc* func, double* obj);

SLEQP_NODISCARD SLEQP_RETCODE
sleqp_func_obj_grad(SleqpFunc* func, SleqpVec* obj_grad);

SLEQP_NODISCARD SLEQP_RETCODE
sleqp_func_cons_val(SleqpFunc* func, SleqpVec* cons_val);

SLEQP_NODISCARD SLEQP_RETCODE
sleqp_func_cons_jac(SleqpFunc* func, SleqpSparseMatrix* cons_jac);

SLEQP_FUNC_FLAGS
sleqp_func_flags(const SleqpFunc* func);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_func_flags_add(SleqpFunc* func, SLEQP_FUNC_FLAGS flags);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_func_flags_remove(SleqpFunc* func, SLEQP_FUNC_FLAGS flags);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_func_flags_set(SleqpFunc* func, SLEQP_FUNC_FLAGS flags, bool value);

bool
sleqp_func_has_flags(const SleqpFunc* func, SLEQP_FUNC_FLAGS flags);

bool
sleqp_func_flags_copy(const SleqpFunc* source,
                      SleqpFunc* target,
                      SLEQP_FUNC_FLAGS flags);

SLEQP_FUNC_TYPE
sleqp_func_get_type(const SleqpFunc* func);

SLEQP_RETCODE
sleqp_func_set_type(SleqpFunc* func, SLEQP_FUNC_TYPE func_type);

/**
 * Returns the setting timer of this function. This timer records
 * the setting of function values.
 *
 * @param[in]     func            The function
 *
 **/
SleqpTimer*
sleqp_func_get_set_timer(SleqpFunc* func);

/**
 * Returns the evaluation timer of this function. This timer records
 * the evaluations of function values.
 *
 * @param[in]     func            The function
 *
 **/
SleqpTimer*
sleqp_func_get_val_timer(SleqpFunc* func);

/**
 * Returns the gradient timer of this function. This timer records
 * the evaluations of function gradients.
 *
 * @param[in]     func            The function
 *
 **/
SleqpTimer*
sleqp_func_get_grad_timer(SleqpFunc* func);

/**
 * Returns the constraint timer of this function. This timer records
 * the evaluations of constraint values.
 *
 * @param[in]     func            The function
 *
 **/
SleqpTimer*
sleqp_func_get_cons_val_timer(SleqpFunc* func);

/**
 * Returns the Jacobian timer of this function. This timer records
 * the evaluations of constraint Jacobians.
 *
 * @param[in]     func            The function
 *
 **/
SleqpTimer*
sleqp_func_get_cons_jac_timer(SleqpFunc* func);

/**
 * Returns the Hessian timer of this function. This timer records
 * the evaluations of Hessian products.
 *
 * @param[in]     func            The function
 *
 **/
SleqpTimer*
sleqp_func_get_hess_timer(SleqpFunc* func);

/**
 * Evaluates the product of the Hessian of the Lagrangian of the given function.
 *
 * @param[in]     func              The function
 * @param[in]     obj_dual          The value \f$ \lambda_0 \f$
 * @param[in]     direction         The direction \f$ d \f$
 * @param[in]     cons_duals        The values \f$ \lambda \f$
 * @param[out]    product           The resulting product
 *
 */
SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_func_hess_prod(SleqpFunc* func,
                     const double* obj_dual,
                     const SleqpVec* direction,
                     const SleqpVec* cons_duals,
                     SleqpVec* product);

/**
 * Evaluates the bilinear product of the Hessian of the Lagrangian of the given
 * function.
 *
 * @param[in]     func              The function
 * @param[in]     obj_dual          The value \f$ \lambda_0 \f$
 * @param[in]     direction         The direction \f$ d \f$
 * @param[in]     cons_duals        The values \f$ \lambda \f$
 * @param[out]    bilinear_prod     The resulting bilinear product
 *
 */
SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_func_hess_bilinear(SleqpFunc* func,
                         const double* obj_dual,
                         const SleqpVec* direction,
                         const SleqpVec* cons_duals,
                         double* bilinear_prod);

#endif /* SLEQP_FUNC_H */
