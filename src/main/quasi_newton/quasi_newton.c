#include "quasi_newton.h"

#include "func.h"
#include "mem.h"

struct SleqpQuasiNewton
{
  int refcount;

  SleqpFunc* func;
  SleqpTimer* update_timer;

  SleqpFunc* quasi_newton_func;

  SleqpQuasiNewtonCallbacks callbacks;
  void* quasi_newton_data;
};

static SLEQP_RETCODE
quasi_newton_func_set_value(SleqpFunc* func,
                            SleqpSparseVec* x,
                            SLEQP_VALUE_REASON reason,
                            bool* reject,
                            int* obj_grad_nnz,
                            int* cons_val_nnz,
                            int* cons_jac_nnz,
                            void* func_data)
{
  SleqpQuasiNewton* quasi_newton = (SleqpQuasiNewton*)func_data;

  SLEQP_CALL(sleqp_func_set_value(quasi_newton->func,
                                  x,
                                  reason,
                                  reject,
                                  obj_grad_nnz,
                                  cons_val_nnz,
                                  cons_jac_nnz));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
quasi_newton_func_obj_val(SleqpFunc* func, double* obj_val, void* func_data)
{
  SleqpQuasiNewton* quasi_newton = (SleqpQuasiNewton*)func_data;

  SLEQP_CALL(sleqp_func_obj_val(quasi_newton->func, obj_val));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
quasi_newton_func_obj_grad(SleqpFunc* func,
                           SleqpSparseVec* obj_grad,
                           void* func_data)
{
  SleqpQuasiNewton* quasi_newton = (SleqpQuasiNewton*)func_data;

  SLEQP_CALL(sleqp_func_obj_grad(quasi_newton->func, obj_grad));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
quasi_newton_func_cons_val(SleqpFunc* func,
                           const SleqpSparseVec* cons_indices,
                           SleqpSparseVec* cons_val,
                           void* func_data)
{
  SleqpQuasiNewton* quasi_newton = (SleqpQuasiNewton*)func_data;

  SLEQP_CALL(sleqp_func_cons_val(quasi_newton->func, cons_indices, cons_val));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
quasi_newton_func_cons_jac(SleqpFunc* func,
                           const SleqpSparseVec* cons_indices,
                           SleqpSparseMatrix* cons_jac,
                           void* func_data)
{
  SleqpQuasiNewton* quasi_newton = (SleqpQuasiNewton*)func_data;

  SLEQP_CALL(sleqp_func_cons_jac(quasi_newton->func, cons_indices, cons_jac));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
quasi_newton_func_hess_prod(SleqpFunc* func,
                            const double* obj_dual,
                            const SleqpSparseVec* direction,
                            const SleqpSparseVec* cons_duals,
                            SleqpSparseVec* product,
                            void* func_data)
{
  SleqpQuasiNewton* quasi_newton = (SleqpQuasiNewton*)func_data;

  SLEQP_CALL(sleqp_quasi_newton_hess_prod(quasi_newton, direction, product));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
quasi_newton_func_create(SleqpQuasiNewton* quasi_newton)
{
  SleqpFunc* func = quasi_newton->func;

  const int num_variables   = sleqp_func_num_vars(func);
  const int num_constraints = sleqp_func_num_cons(func);

  SleqpFuncCallbacks callbacks = {.set_value = quasi_newton_func_set_value,
                                  .obj_val   = quasi_newton_func_obj_val,
                                  .obj_grad  = quasi_newton_func_obj_grad,
                                  .cons_val  = quasi_newton_func_cons_val,
                                  .cons_jac  = quasi_newton_func_cons_jac,
                                  .hess_prod = quasi_newton_func_hess_prod,
                                  .func_free = NULL};

  SLEQP_CALL(sleqp_func_create(&quasi_newton->quasi_newton_func,
                               &callbacks,
                               num_variables,
                               num_constraints,
                               quasi_newton));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_quasi_newton_create(SleqpQuasiNewton** star,
                          SleqpFunc* func,
                          SleqpQuasiNewtonCallbacks* callbacks,
                          void* quasi_newton_data)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpQuasiNewton* quasi_newton = *star;

  *quasi_newton = (SleqpQuasiNewton){0};

  quasi_newton->refcount = 1;

  SLEQP_CALL(sleqp_func_capture(func));
  quasi_newton->func = func;

  SLEQP_CALL(sleqp_timer_create(&(quasi_newton->update_timer)));

  SLEQP_CALL(quasi_newton_func_create(quasi_newton));

  quasi_newton->callbacks         = *callbacks;
  quasi_newton->quasi_newton_data = quasi_newton_data;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_quasi_newton_push(SleqpQuasiNewton* quasi_newton,
                        const SleqpIterate* old_iterate,
                        const SleqpIterate* new_iterate,
                        const SleqpSparseVec* multipliers)
{
  SLEQP_CALL(sleqp_timer_start(quasi_newton->update_timer));

  SLEQP_CALL(quasi_newton->callbacks.push(old_iterate,
                                          new_iterate,
                                          multipliers,
                                          quasi_newton->quasi_newton_data));

  SLEQP_CALL(sleqp_timer_stop(quasi_newton->update_timer));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_quasi_newton_reset(SleqpQuasiNewton* quasi_newton)
{
  SLEQP_CALL(quasi_newton->callbacks.reset(quasi_newton->quasi_newton_data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_quasi_newton_hess_prod(SleqpQuasiNewton* quasi_newton,
                             const SleqpSparseVec* direction,
                             SleqpSparseVec* product)
{
  SLEQP_CALL(quasi_newton->callbacks
               .hess_prod(direction, product, quasi_newton->quasi_newton_data));

  return SLEQP_OKAY;
}

SleqpTimer*
sleqp_quasi_newton_update_timer(SleqpQuasiNewton* quasi_newton)
{
  return quasi_newton->update_timer;
}

SleqpFunc*
sleqp_quasi_newton_get_func(SleqpQuasiNewton* quasi_newton)
{
  return quasi_newton->quasi_newton_func;
}

SLEQP_RETCODE
sleqp_quasi_newton_capture(SleqpQuasiNewton* quasi_newton)
{
  ++quasi_newton->refcount;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
quasi_newton_free(SleqpQuasiNewton** star)
{
  SleqpQuasiNewton* quasi_newton = *star;

  SLEQP_CALL(sleqp_timer_free(&(quasi_newton->update_timer)));

  SLEQP_CALL(quasi_newton->callbacks.free(quasi_newton->quasi_newton_data));

  SLEQP_CALL(sleqp_func_release(&(quasi_newton->func)));

  SLEQP_CALL(sleqp_func_release(&quasi_newton->quasi_newton_func));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_quasi_newton_release(SleqpQuasiNewton** star)
{
  SleqpQuasiNewton* quasi_newton = *star;

  if (!quasi_newton)
  {
    return SLEQP_OKAY;
  }

  if (--(quasi_newton->refcount) == 0)
  {
    SLEQP_CALL(quasi_newton_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_quasi_newton_create_default(SleqpQuasiNewton** star,
                                  SleqpFunc* func,
                                  SleqpParams* params,
                                  SleqpOptions* options)
{
  const SLEQP_HESS_EVAL hessian_eval
    = sleqp_options_int_value(options, SLEQP_OPTION_INT_HESS_EVAL);

  if (hessian_eval == SLEQP_HESS_EVAL_SIMPLE_BFGS
      || hessian_eval == SLEQP_HESS_EVAL_DAMPED_BFGS)
  {
    SLEQP_CALL(sleqp_bfgs_create(star, func, params, options));
  }
  else if (hessian_eval == SLEQP_HESS_EVAL_SR1)
  {
    SLEQP_CALL(sleqp_sr1_create(star, func, params, options));
  }
  else
  {
    *star = NULL;
  }

  return SLEQP_OKAY;
}
