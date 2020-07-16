#include "sleqp_func.h"

#include "sleqp_mem.h"

struct SleqpFunc
{
  int refcount;

  SleqpFuncCallbacks callbacks;

  int num_variables;
  void* data;

  int num_func_evals;
  int num_cons_evals;

  int num_grad_evals;
  int num_jac_evals;

  int num_hess_evals;

  SleqpTimer* eval_timer;
  SleqpTimer* hess_timer;

  SleqpSparseVec* product;

  SleqpHessianStruct* hess_struct;
};

SLEQP_RETCODE sleqp_func_create(SleqpFunc** fstar,
                                SleqpFuncCallbacks* callbacks,
                                int num_variables,
                                void* func_data)
{
  SLEQP_CALL(sleqp_malloc(fstar));

  SleqpFunc* func = *fstar;

  *func = (SleqpFunc) {0};
  func->refcount = 1;

  func->callbacks = *callbacks;

  func->num_variables = num_variables;
  func->data = func_data;

  SLEQP_CALL(sleqp_timer_create(&func->eval_timer));
  SLEQP_CALL(sleqp_timer_create(&func->hess_timer));

  SLEQP_CALL(sleqp_sparse_vector_create(&func->product, num_variables, 0));

  SLEQP_CALL(sleqp_hessian_struct_create(&func->hess_struct,
                                         num_variables,
                                         false));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_func_set_value(SleqpFunc* func,
                                   SleqpSparseVec* x,
                                   SLEQP_VALUE_REASON reason,
                                   int* func_grad_nnz,
                                   int* cons_val_nnz,
                                   int* cons_jac_nnz)
{
  SLEQP_CALL(func->callbacks.set_value(x,
                                       reason,
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

  SLEQP_CALL(sleqp_timer_start(func->eval_timer));

  SLEQP_CALL(func->callbacks.func_eval(func->num_variables,
                                       cons_indices,
                                       func_val,
                                       func_grad,
                                       cons_val,
                                       cons_jac,
                                       func->data));

  SLEQP_CALL(sleqp_timer_stop(func->eval_timer));

  return SLEQP_OKAY;
}

SleqpHessianStruct* sleqp_func_get_hess_struct(SleqpFunc* func)
{
  return func->hess_struct;
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

SleqpTimer* sleqp_func_get_eval_timer(SleqpFunc* func)
{
  return func->eval_timer;
}

SleqpTimer* sleqp_func_get_hess_timer(SleqpFunc* func)
{
  return func->hess_timer;
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

  SLEQP_CALL(sleqp_timer_start(func->hess_timer));

  SLEQP_CALL(func->callbacks.hess_prod(func->num_variables,
                                       func_dual,
                                       direction,
                                       cons_duals,
                                       product,
                                       func->data));

  SLEQP_CALL(sleqp_timer_stop(func->hess_timer));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_func_hess_bilinear(SleqpFunc* func,
                                       double* func_dual,
                                       SleqpSparseVec* direction,
                                       SleqpSparseVec* cons_duals,
                                       double* bilinear_prod)
{
  SLEQP_CALL(func->callbacks.hess_prod(func->num_variables,
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

static SLEQP_RETCODE func_free(SleqpFunc** fstar)
{
  SleqpFunc* func = *fstar;

  if(!func)
  {
    return SLEQP_OKAY;
  }

  if(func->callbacks.func_free)
  {
    SLEQP_CALL(func->callbacks.func_free(func->data));
  }

  SLEQP_CALL(sleqp_timer_free(&func->hess_timer));
  SLEQP_CALL(sleqp_timer_free(&func->eval_timer));

  SLEQP_CALL(sleqp_sparse_vector_free(&func->product));

  SLEQP_CALL(sleqp_hessian_struct_free(&func->hess_struct));

  sleqp_free(fstar);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_func_capture(SleqpFunc* func)
{
  ++func->refcount;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_func_release(SleqpFunc** fstar)
{
  SleqpFunc* func = *fstar;

  if(!func)
  {
    return SLEQP_OKAY;
  }

  if(--func->refcount == 0)
  {
    SLEQP_CALL(func_free(fstar));
  }

  *fstar = NULL;

  return SLEQP_OKAY;
}
