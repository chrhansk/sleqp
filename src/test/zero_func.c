#include "zero_func.h"

#include "lsq.h"

static SLEQP_RETCODE
zero_func_set(SleqpFunc* func,
              SleqpSparseVec* x,
              SLEQP_VALUE_REASON reason,
              bool* reject,
              int* obj_grad_nnz,
              int* cons_val_nnz,
              int* cons_jac_nnz,
              void* func_data)
{
  *obj_grad_nnz = 0;
  *cons_val_nnz = 0;
  *cons_jac_nnz = 0;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
zero_func_obj_val(SleqpFunc* func, double* obj_val, void* func_data)
{
  *obj_val = 0.;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
zero_func_obj_grad(SleqpFunc* func, SleqpSparseVec* obj_grad, void* func_data)
{
  return SLEQP_OKAY;
}

static SLEQP_RETCODE
zero_func_cons_val(SleqpFunc* func, SleqpSparseVec* cons_val, void* func_data)
{
  return SLEQP_OKAY;
}

static SLEQP_RETCODE
zero_func_cons_jac(SleqpFunc* func,
                   SleqpSparseMatrix* cons_jac,
                   void* func_data)
{
  return SLEQP_OKAY;
}

static SLEQP_RETCODE
zero_func_hess_prod(SleqpFunc* func,
                    const double* obj_dual,
                    const SleqpSparseVec* direction,
                    const SleqpSparseVec* cons_duals,
                    SleqpSparseVec* result,
                    void* func_data)
{
  return SLEQP_OKAY;
}

static SLEQP_RETCODE
zero_lsq_func_residuals(SleqpFunc* func,
                        SleqpSparseVec* residual,
                        void* func_data)
{
  return SLEQP_OKAY;
}

static SLEQP_RETCODE
zero_lsq_func_jac_forward(SleqpFunc* func,
                          const SleqpSparseVec* forward_direction,
                          SleqpSparseVec* product,
                          void* func_data)
{
  return SLEQP_OKAY;
}

static SLEQP_RETCODE
zero_lsq_func_jac_adjoint(SleqpFunc* func,
                          const SleqpSparseVec* adjoint_direction,
                          SleqpSparseVec* product,
                          void* func_data)
{
  return SLEQP_OKAY;
}

SLEQP_RETCODE
zero_func_create(SleqpFunc** star, int num_variables, int num_constraints)
{

  SleqpFuncCallbacks callbacks = {.set_value = zero_func_set,
                                  .obj_val   = zero_func_obj_val,
                                  .obj_grad  = zero_func_obj_grad,
                                  .cons_val  = zero_func_cons_val,
                                  .cons_jac  = zero_func_cons_jac,
                                  .hess_prod = zero_func_hess_prod,
                                  .func_free = NULL};

  SLEQP_CALL(
    sleqp_func_create(star, &callbacks, num_variables, num_constraints, NULL));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
zero_lsq_func_create(SleqpFunc** star,
                     SleqpParams* params,
                     int num_variables,
                     int num_constraints,
                     int num_residuals)
{
  SleqpLSQCallbacks callbacks = {.set_value       = zero_func_set,
                                 .lsq_residuals   = zero_lsq_func_residuals,
                                 .lsq_jac_forward = zero_lsq_func_jac_forward,
                                 .lsq_jac_adjoint = zero_lsq_func_jac_adjoint,
                                 .cons_val        = zero_func_cons_val,
                                 .cons_jac        = zero_func_cons_jac,
                                 .func_free       = NULL};

  SLEQP_CALL(sleqp_lsq_func_create(star,
                                   &callbacks,
                                   num_variables,
                                   num_constraints,
                                   num_residuals,
                                   0.,
                                   params,
                                   NULL));

  return SLEQP_OKAY;
}
