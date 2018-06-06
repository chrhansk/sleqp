#include "sleqp_cauchy.h"

SLEQP_RETCODE sleqp_cauchy_direction(SleqpProblem* problem,
                                     SleqpSparseVec* x,
                                     SleqpLPi* lp_interface,
                                     SleqpSparseVec* func_grad,
                                     SleqpSparseMatrix* cons_jac)
{
  double func_val;

  SLEQP_CALL(sleqp_func_eval(problem->func,
                             problem->num_variables,
                             NULL,
                             &func_val,
                             func_grad,
                             cons_jac));

  return SLEQP_OKAY;
}
