#include "sleqp_func.h"

#include "sleqp_mem.h"

struct SleqpFunc
{
  SLEQP_FUNC_SET set_value;
  SLEQP_FUNC_EVAL eval;
  SLEQP_HESS_PRODUCT eval_hess_prod;

  int num_variables;
  void* data;

  SleqpSparseVec* product;
};

SLEQP_RETCODE sleqp_func_create(SleqpFunc** fstar,
                                SLEQP_FUNC_SET set_value,
                                SLEQP_FUNC_EVAL eval,
                                SLEQP_HESS_PRODUCT eval_hess_prod,
                                int num_variables,
                                void* func_data)
{
  SLEQP_CALL(sleqp_malloc(fstar));

  SleqpFunc* func = *fstar;

  func->set_value = set_value;
  func->eval = eval;
  func->eval_hess_prod = eval_hess_prod;

  func->num_variables = num_variables;
  func->data = func_data;

  SLEQP_CALL(sleqp_sparse_vector_create(&func->product, num_variables, 0));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_func_set_value(SleqpFunc* func,
                                   SleqpSparseVec* x,
                                   int* func_grad_nnz,
                                   int* cons_val_nnz,
                                   int* cons_jac_nnz)
{


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


  SLEQP_CALL(func->eval(func->num_variables,
                        indices,
                        func_val,
                        func_grad,
                        cons_val,
                        cons_jac,
                        func->data));

  return SLEQP_OKAY;
}

int sleqp_func_get_num_variables(SleqpFunc* func)
{
  return func->num_variables;
}

SLEQP_RETCODE sleqp_func_hess_product(SleqpFunc* func,
                                      double* func_dual,
                                      SleqpSparseVec* direction,
                                      SleqpSparseVec* cons_duals,
                                      SleqpSparseVec* product)
{
  SLEQP_CALL(func->eval_hess_prod(func->num_variables,
                                  func_dual,
                                  direction,
                                  cons_duals,
                                  product,
                                  func->data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_func_hess_bilinear(SleqpFunc* func,
                                       double* func_dual,
                                       SleqpSparseVec* direction,
                                       SleqpSparseVec* cons_duals,
                                       double* bilinear_prod)
{
  SLEQP_CALL(func->eval_hess_prod(func->num_variables,
                                  func_dual,
                                  direction,
                                  cons_duals,
                                  func->product,
                                  func->data));

  SLEQP_CALL(sleqp_sparse_vector_dot(direction,
                                     func->product,
                                     bilinear_prod));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_func_free(SleqpFunc** fstar)
{
  SleqpFunc* func = *fstar;

  SLEQP_CALL(sleqp_sparse_vector_free(&func->product));

  sleqp_free(fstar);

  return SLEQP_OKAY;
}
