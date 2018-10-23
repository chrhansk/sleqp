#include "sleqp.h"

#include <assert.h>

SLEQP_RETCODE sleqp_set_and_evaluate(SleqpProblem* problem,
                                     SleqpIterate* iterate)
{
  int func_grad_nnz = 0;
  int cons_val_nnz = 0;
  int cons_jac_nnz = 0;

  SLEQP_CALL(sleqp_func_set_value(problem->func,
                                  iterate->x,
                                  &func_grad_nnz,
                                  &cons_val_nnz,
                                  &cons_jac_nnz));

  SLEQP_CALL(sleqp_sparse_vector_reserve(iterate->func_grad, func_grad_nnz));

  SLEQP_CALL(sleqp_sparse_vector_reserve(iterate->cons_val, cons_val_nnz));

  SLEQP_CALL(sleqp_sparse_matrix_reserve(iterate->cons_jac, cons_jac_nnz));

  SLEQP_CALL(sleqp_func_eval(problem->func,
                             NULL,
                             &iterate->func_val,
                             iterate->func_grad,
                             iterate->cons_val,
                             iterate->cons_jac));

  assert(sleqp_sparse_vector_valid(iterate->func_grad));
  assert(sleqp_sparse_vector_valid(iterate->cons_val));

  assert(sleqp_sparse_matrix_valid(iterate->cons_jac));

  return SLEQP_OKAY;
}
