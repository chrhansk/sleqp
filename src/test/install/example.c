#include <stdlib.h>

#include <sleqp.h>

#define MAIN_CALL(x)                            \
  do                                            \
  {                                             \
    SLEQP_RETCODE _retcode_ = (x);              \
    if(_retcode_ != SLEQP_OKAY)                 \
    {                                           \
      return EXIT_FAILURE;                      \
    }                                           \
  }                                             \
  while(0)

SLEQP_RETCODE test_set(SleqpFunc* func,
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

SLEQP_RETCODE test_val(SleqpFunc* func,
                             double* func_val,
                             void* func_data)
{
  *func_val = 0.;

  return SLEQP_OKAY;
}

SLEQP_RETCODE test_grad(SleqpFunc* func,
                        SleqpSparseVec* func_grad,
                        void* func_data)
{
  return SLEQP_OKAY;
}

SLEQP_RETCODE test_hess_prod(SleqpFunc* func,
                                   const double* func_dual,
                                   const SleqpSparseVec* direction,
                                   const SleqpSparseVec* cons_duals,
                                   SleqpSparseVec* product,
                                   void* func_data)
{
  return SLEQP_OKAY;
}

int main(int argc, char *argv[])
{
  SleqpFunc* func;

  const int num_variables = 1;
  const int num_constraints = 0;

  SleqpFuncCallbacks callbacks = {
    .set_value = test_set,
    .func_val  = test_val,
    .func_grad = test_grad,
    .cons_val  = NULL,
    .cons_jac  = NULL,
    .hess_prod = test_hess_prod,
    .func_free = NULL
  };

  MAIN_CALL(sleqp_func_create(&func,
                              &callbacks,
                              num_variables,
                              num_constraints,
                              NULL));

  MAIN_CALL(sleqp_func_release(&func));

  return EXIT_SUCCESS;
}
