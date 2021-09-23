#include "zero_func.h"

static SLEQP_RETCODE
zero_func_set(SleqpFunc* func,
              SleqpSparseVec* x,
              SLEQP_VALUE_REASON reason,
              bool* reject,
              int* func_grad_nnz,
              int* cons_val_nnz,
              int* cons_jac_nnz,
              void* func_data)
{
  *func_grad_nnz = 0;
  *cons_val_nnz = 0;
  *cons_jac_nnz = 0;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
zero_func_val(SleqpFunc* func,
              double* func_val,
              void* func_data)
{
  *func_val = 0.;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
zero_func_grad(SleqpFunc* func,
               SleqpSparseVec* func_grad,
               void* func_data)
{
  return SLEQP_OKAY;
}

static SLEQP_RETCODE
zero_func_cons_val(SleqpFunc* func,
                   const SleqpSparseVec* cons_indices,
                   SleqpSparseVec* cons_val,
                   void* func_data)
{
  return SLEQP_OKAY;
}

static SLEQP_RETCODE
zero_func_cons_jac(SleqpFunc* func,
                   const SleqpSparseVec* cons_indices,
                   SleqpSparseMatrix* cons_jac,
                   void* func_data)
{
  return SLEQP_OKAY;
}

static SLEQP_RETCODE
zero_func_hess_prod(SleqpFunc* func,
                    const double* func_dual,
                    const SleqpSparseVec* direction,
                    const SleqpSparseVec* cons_duals,
                    SleqpSparseVec* result,
                    void* func_data)
{
  return SLEQP_OKAY;
}




SLEQP_RETCODE zero_func_create(SleqpFunc** star,
                               int num_variables,
                               int num_constraints)
{

  SleqpFuncCallbacks callbacks = {
    .set_value = zero_func_set,
    .func_val  = zero_func_val,
    .func_grad = zero_func_grad,
    .cons_val  = zero_func_cons_val,
    .cons_jac  = zero_func_cons_jac,
    .hess_prod = zero_func_hess_prod,
    .func_free = NULL
  };

  SLEQP_CALL(sleqp_func_create(star,
                               &callbacks,
                               num_variables,
                               num_constraints,
                               NULL));

  return SLEQP_OKAY;
}
