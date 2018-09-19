#include "sleqp_func.h"

#include "sleqp_mem.h"

struct SleqpFunc
{
  SLEQP_FUNC_SET set_value;
  SLEQP_FUNC_EVAL eval;
  SLEQP_HESS_EVAL_BILINEAR eval_bilin;

  size_t num_variables;
  void* data;
};

SLEQP_RETCODE sleqp_func_create(SleqpFunc** fstar,
                                SLEQP_FUNC_SET set_value,
                                SLEQP_FUNC_EVAL eval,
                                SLEQP_HESS_EVAL_BILINEAR eval_bilin,
                                size_t num_variables,
                                void* func_data)
{
  sleqp_malloc(fstar);

  SleqpFunc* func = *fstar;

  func->set_value = set_value;
  func->eval = eval;
  func->eval_bilin = eval_bilin;

  func->num_variables = num_variables;
  func->data = func_data;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_func_set_value(SleqpFunc* func,
                                   SleqpSparseVec* x,
                                   size_t* func_grad_nnz,
                                   size_t* cons_val_nnz,
                                   size_t* cons_jac_nnz)
{
  assert(func);

  SLEQP_CALL(func->set_value(x,
                             func->num_variables,
                             func_grad_nnz,
                             cons_val_nnz,
                             cons_jac_nnz,
                             func->data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_func_eval(SleqpFunc* func,
                              int* indices,
                              double* func_val,
                              SleqpSparseVec* func_grad,
                              SleqpSparseVec* cons_val,
                              SleqpSparseMatrix* cons_jac)
{
  assert(func);

  SLEQP_CALL(func->eval(func->num_variables,
                        indices,
                        func_val,
                        func_grad,
                        cons_val,
                        cons_jac,
                        func->data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_hess_eval_bilinear(SleqpFunc* func,
                                       double* func_dual,
                                       SleqpSparseVec* direction,
                                       SleqpSparseVec* cons_duals,
                                       double* bilinear_prod,
                                       void* func_data)
{
  assert(func);

  SLEQP_CALL(func->eval_bilin(func->num_variables,
                              func_dual,
                              direction,
                              cons_duals,
                              bilinear_prod,
                              func->data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_func_free(SleqpFunc** fstar)
{
  sleqp_free(fstar);

  return SLEQP_OKAY;
}
