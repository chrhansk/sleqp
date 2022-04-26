#ifndef SLEQP_PUB_DYN_H
#define SLEQP_PUB_DYN_H

#include "sleqp/pub_func.h"

/**
 * Evaluates the dynamic function at the current input vector
 *
 * @param[in]     func            The function
 * @param[in]     accuracy        The desired accuracy \f$  \epsilon \f$
 * @param[out]    obj_val         The function value \f$ f(x, \epsilon) \f$
 * @param[in,out] func_data       The function data
 **/
typedef SLEQP_RETCODE (*SLEQP_DYN_FUNC_OBJ_VAL)(SleqpFunc* func,
                                                double accuracy,
                                                double* obj_val,
                                                void* func_data);

/**
 * Evaluates the dynamic constraints at the current input vector
 *
 * @param[in]     func            The function
 * @param[in]     accuracy        The desired accuracy \f$  \epsilon \f$
 * @param[out]    cons_val        The value of the constraint function \f$ c(x,
 *\epsilon) \f$
 * @param[in,out] func_data       The function data
 **/
typedef SLEQP_RETCODE (*SLEQP_DYN_FUNC_CONS_VAL)(SleqpFunc* func,
                                                 double accuracy,
                                                 SleqpVec* cons_val,
                                                 void* func_data);

typedef struct
{
  SLEQP_FUNC_SET set_value;
  SLEQP_FUNC_NONZEROS nonzeros;
  SLEQP_DYN_FUNC_OBJ_VAL obj_val;
  SLEQP_FUNC_OBJ_GRAD obj_grad;
  SLEQP_DYN_FUNC_CONS_VAL cons_val;
  SLEQP_FUNC_CONS_JAC cons_jac;
  SLEQP_FUNC_HESS_PROD hess_prod;
  SLEQP_FUNC_FREE func_free;
} SleqpDynFuncCallbacks;

/**
 * Creates a new function
 *
 * @param[out] fstar            A pointer to the function to be created
 * @param[in]  callbacks        A callback to the function callbacks
 * @param[in]  num_variables    The number of variables
 * @param[in]  num_constraints  The number of constraints
 * @param[in]  func_data        The function data
 **/
SLEQP_EXPORT SLEQP_NODISCARD SLEQP_RETCODE
sleqp_dyn_func_create(SleqpFunc** fstar,
                      SleqpDynFuncCallbacks* callbacks,
                      int num_variables,
                      int num_constraints,
                      void* func_data);

#endif /* SLEQP_PUB_DYN_H */
