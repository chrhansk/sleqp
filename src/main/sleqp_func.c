#include "sleqp_func.h"

#include "sleqp_mem.h"

struct SleqpFunc
{
  SLEQP_FUNC_SET set_value;
  SLEQP_FUNC_EVAL eval;
  SLEQP_HESS_PRODUCT eval_hess_prod;

  int num_variables;
  void* data;

  int num_func_evals;
  int num_cons_evals;

  int num_grad_evals;
  int num_jac_evals;

  int num_hess_evals;

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

  func->num_func_evals = 0;
  func->num_cons_evals = 0;

  func->num_grad_evals = 0;
  func->num_jac_evals = 0;

  func->num_hess_evals = 0;

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
                              SleqpSparseVec* cons_indices,
                              double* func_val,
                              SleqpSparseVec* func_grad,
                              SleqpSparseVec* cons_val,
                              SleqpSparseMatrix* cons_jac)
{
  if(func_grad)
  {
    ++func->num_grad_evals;

    SLEQP_CALL(sleqp_sparse_vector_clear(func_grad));
  }

  if(cons_val)
  {
    ++func->num_cons_evals;

    SLEQP_CALL(sleqp_sparse_vector_clear(cons_val));
  }

  if(cons_jac)
  {
    ++func->num_jac_evals;

    SLEQP_CALL(sleqp_sparse_matrix_clear(cons_jac));
  }


  ++func->num_func_evals;

  SLEQP_CALL(func->eval(func->num_variables,
                        cons_indices,
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

int sleqp_func_get_num_func_evals(SleqpFunc* func)
{
  return func->num_func_evals;
}
int sleqp_func_get_num_cons_evals(SleqpFunc* func)
{
  return func->num_cons_evals;
}

int sleqp_func_get_num_grad_evals(SleqpFunc* func)
{
  return func->num_grad_evals;
}
int sleqp_func_get_num_jac_evals(SleqpFunc* func)
{
  return func->num_jac_evals;
}

int sleqp_func_get_num_hess_evals(SleqpFunc* func)
{
  return func->num_hess_evals;
}

SLEQP_RETCODE sleqp_func_hess_prod(SleqpFunc* func,
                                   double* func_dual,
                                   SleqpSparseVec* direction,
                                   SleqpSparseVec* cons_duals,
                                   SleqpSparseVec* product)
{
  assert(func->num_variables == direction->dim);
  assert(func->num_variables == product->dim);

  ++func->num_hess_evals;

  SLEQP_CALL(sleqp_sparse_vector_clear(product));

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

void* sleqp_func_get_data(SleqpFunc* func)
{
  return func->data;
}

SLEQP_RETCODE sleqp_func_free(SleqpFunc** fstar)
{
  SleqpFunc* func = *fstar;

  if(!func)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_sparse_vector_free(&func->product));

  sleqp_free(fstar);

  return SLEQP_OKAY;
}
