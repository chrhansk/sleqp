#include <stdlib.h>

#include <sleqp.h>

#define MAIN_CALL(x)                                                           \
  do                                                                           \
  {                                                                            \
    SLEQP_RETCODE _retcode_ = (x);                                             \
    if (_retcode_ != SLEQP_OKAY)                                               \
    {                                                                          \
      return EXIT_FAILURE;                                                     \
    }                                                                          \
  } while (0)

SLEQP_RETCODE
test_set(SleqpFunc* func,
         SleqpVec* x,
         SLEQP_VALUE_REASON reason,
         bool* reject,
         void* func_data)
{
  return SLEQP_OKAY;
}

SLEQP_RETCODE
test_obj_val(SleqpFunc* func, double* obj, void* func_data)
{
  *obj = 0.;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
test_obj_grad(SleqpFunc* func, SleqpVec* obj_grad, void* func_data)
{
  return SLEQP_OKAY;
}

SLEQP_RETCODE
test_hess_prod(SleqpFunc* func,
               const double* obj_dual,
               const SleqpVec* direction,
               const SleqpVec* cons_duals,
               SleqpVec* product,
               void* func_data)
{
  return SLEQP_OKAY;
}

int
main(int argc, char* argv[])
{
  SleqpFunc* func;

  const int num_variables   = 1;
  const int num_constraints = 0;

  SleqpFuncCallbacks callbacks = {.set_value = test_set,
                                  .obj_val   = test_obj_val,
                                  .obj_grad  = test_obj_grad,
                                  .cons_val  = NULL,
                                  .cons_jac  = NULL,
                                  .hess_prod = test_hess_prod,
                                  .func_free = NULL};

  MAIN_CALL(
    sleqp_func_create(&func, &callbacks, num_variables, num_constraints, NULL));

  MAIN_CALL(sleqp_func_release(&func));

  return EXIT_SUCCESS;
}
