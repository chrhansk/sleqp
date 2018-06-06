#include "sleqp_func.h"

#include "sleqp_mem.h"

struct SleqpFunc
{
  SLEQP_FUNC_SET set_value;
  SLEQP_FUNC_EVAL eval;
  SLEQP_HESS_EVAL_BILINEAR eval_bilin;
  void* data;
};

SLEQP_RETCODE sleqp_func_create(SleqpFunc** fstar,
                                SLEQP_FUNC_SET set_value,
                                SLEQP_FUNC_EVAL eval,
                                SLEQP_HESS_EVAL_BILINEAR eval_bilin,
                                void* func_data)
{
  sleqp_malloc(fstar);

  SleqpFunc* func = *fstar;

  func->set_value = set_value;
  func->eval = eval;
  func->eval_bilin = eval_bilin;
  func->data = func_data;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_func_set_value(SleqpFunc* func,
                                   SleqpSparseVec* x,
                                   size_t num_variables,
                                   size_t* grad_nnz,
                                   size_t* jac_nnz)
{
  assert(func);

  SLEQP_CALL(func->set_value(x,
                             num_variables,
                             grad_nnz,
                             jac_nnz,
                             func->data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_func_eval(SleqpFunc* func,
                              size_t num_variables,
                              int* indices,
                              double* func_val,
                              SleqpSparseVec* func_grad,
                              SleqpSparseMatrix* cons_jac)
{
  assert(func);

  SLEQP_CALL(func->eval(num_variables,
                        indices,
                        func_val,
                        func_grad,
                        cons_jac,
                        func->data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_hess_eval_bilinear(SleqpFunc* func,
                                       size_t num_variables,
                                       double* fval,
                                       SleqpSparseVec* direction,
                                       SleqpSparseVec* multipliers)
{
  assert(func);

  SLEQP_CALL(func->eval_bilin(num_variables,
                              fval,
                              direction,
                              multipliers,
                              func->data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_func_free(SleqpFunc** fstar)
{
  sleqp_free(fstar);

  return SLEQP_OKAY;
}
