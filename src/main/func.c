#include "func.h"

#include "cmp.h"
#include "fail.h"
#include "log.h"
#include "mem.h"

#include "sparse/sparse_matrix.h"

struct SleqpFunc
{
  int refcount;

  SleqpFuncCallbacks callbacks;

  bool hessian_psd;
  SLEQP_FUNC_TYPE type;

  int num_variables;
  int num_constraints;

  void* data;

  SleqpTimer* set_timer;
  SleqpTimer* val_timer;
  SleqpTimer* grad_timer;

  SleqpTimer* cons_val_timer;
  SleqpTimer* cons_jac_timer;

  SleqpTimer* hess_timer;

  SleqpSparseVec* product;

  SleqpHessianStruct* hess_struct;
};

SLEQP_RETCODE sleqp_func_create(SleqpFunc** fstar,
                                SleqpFuncCallbacks* callbacks,
                                int num_variables,
                                int num_constraints,
                                void* func_data)
{
  SLEQP_CALL(sleqp_malloc(fstar));

  SleqpFunc* func = *fstar;

  *func = (SleqpFunc) {0};
  func->refcount = 1;

  func->callbacks = *callbacks;
  func->hessian_psd = false;

  func->num_variables = num_variables;
  func->num_constraints = num_constraints;
  func->data = func_data;
  func->type = SLEQP_FUNC_TYPE_REGULAR;

  SLEQP_CALL(sleqp_timer_create(&func->set_timer));
  SLEQP_CALL(sleqp_timer_create(&func->val_timer));
  SLEQP_CALL(sleqp_timer_create(&func->grad_timer));

  SLEQP_CALL(sleqp_timer_create(&func->cons_val_timer));
  SLEQP_CALL(sleqp_timer_create(&func->cons_jac_timer));

  SLEQP_CALL(sleqp_timer_create(&func->hess_timer));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&func->product, num_variables));

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
  assert(sleqp_sparse_vector_is_valid(x));
  assert(sleqp_sparse_vector_is_finite(x));

  SLEQP_CALL(sleqp_timer_start(func->set_timer));

  SLEQP_CALL(func->callbacks.set_value(func,
                                       x,
                                       reason,
                                       func_grad_nnz,
                                       cons_val_nnz,
                                       cons_jac_nnz,
                                       func->data));

  SLEQP_CALL(sleqp_timer_stop(func->set_timer));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_func_val(SleqpFunc* func,
                             double* func_val)
{
  if(func_val)
  {
    SLEQP_CALL(sleqp_timer_start(func->val_timer));

    SLEQP_CALL(func->callbacks.func_val(func,
                                        func_val,
                                        func->data));

    SLEQP_CALL(sleqp_timer_stop(func->val_timer));

    sleqp_assert_msg(sleqp_is_finite(*func_val),
                     "Returned infinite function value");
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_func_grad(SleqpFunc* func,
                              SleqpSparseVec* func_grad)
{
  if(func_grad)
  {
    SLEQP_CALL(sleqp_sparse_vector_clear(func_grad));

    SLEQP_CALL(sleqp_timer_start(func->grad_timer));

    SLEQP_CALL(func->callbacks.func_grad(func,
                                         func_grad,
                                         func->data));

    SLEQP_CALL(sleqp_timer_stop(func->grad_timer));

    sleqp_assert_msg(sleqp_sparse_vector_is_valid(func_grad),
                     "Returned invalid function gradient");

    sleqp_assert_msg(sleqp_sparse_vector_is_finite(func_grad),
                     "Returned function gradient is not all-finite");
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_func_cons_val(SleqpFunc* func,
                                  const SleqpSparseVec* cons_indices,
                                  SleqpSparseVec* cons_val)
{
  const int num_constraints = sleqp_func_get_num_constraints(func);

  if(cons_val)
  {
    SLEQP_CALL(sleqp_sparse_vector_clear(cons_val));

    if((num_constraints != 0) && (func->callbacks.cons_val))
    {
      SLEQP_CALL(sleqp_timer_start(func->cons_val_timer));

      SLEQP_CALL(func->callbacks.cons_val(func,
                                          cons_indices,
                                          cons_val,
                                          func->data));

      SLEQP_CALL(sleqp_timer_stop(func->cons_val_timer));
    }

    sleqp_assert_msg(sleqp_sparse_vector_is_valid(cons_val),
                     "Returned invalid constraint values");

    sleqp_assert_msg(sleqp_sparse_vector_is_finite(cons_val),
                     "Returned constraint values are not all-finite");
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_func_cons_jac(SleqpFunc* func,
                                  const SleqpSparseVec* cons_indices,
                                  SleqpSparseMatrix* cons_jac)
{
  const int num_constraints = sleqp_func_get_num_constraints(func);

  if(cons_jac)
  {
    SLEQP_CALL(sleqp_sparse_matrix_clear(cons_jac));

    if((num_constraints != 0) && (func->callbacks.cons_jac))
    {
      SLEQP_CALL(sleqp_timer_start(func->cons_jac_timer));

      SLEQP_CALL(func->callbacks.cons_jac(func,
                                          cons_indices,
                                          cons_jac,
                                          func->data));

      SLEQP_CALL(sleqp_timer_stop(func->cons_jac_timer));
    }

    sleqp_assert_msg(sleqp_sparse_matrix_is_valid(cons_jac),
                     "Returned invalid constraint Jacobian");

    sleqp_assert_msg(sleqp_sparse_matrix_is_finite(cons_jac),
                     "Returned constraint Jacobian is not all-finite");
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_func_eval(SleqpFunc* func,
                              const SleqpSparseVec* cons_indices,
                              double* func_val,
                              SleqpSparseVec* func_grad,
                              SleqpSparseVec* cons_val,
                              SleqpSparseMatrix* cons_jac)
{
  SLEQP_CALL(sleqp_func_val(func, func_val));

  SLEQP_CALL(sleqp_func_grad(func, func_grad));

  SLEQP_CALL(sleqp_func_cons_val(func, cons_indices, cons_val));

  SLEQP_CALL(sleqp_func_cons_jac(func, cons_indices, cons_jac));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_func_set_callbacks(SleqpFunc* func,
                                       SleqpFuncCallbacks* callbacks)
{
  func->callbacks = *callbacks;

  return SLEQP_OKAY;
}

SleqpHessianStruct* sleqp_func_get_hess_struct(SleqpFunc* func)
{
  return func->hess_struct;
}

bool sleqp_func_has_psd_hessian(SleqpFunc* func)
{
  return func->hessian_psd;
}

SLEQP_RETCODE sleqp_func_set_psd_hessian(SleqpFunc* func,
                                         bool value)
{
  func->hessian_psd = value;

  return SLEQP_OKAY;
}

SLEQP_FUNC_TYPE sleqp_func_get_type(SleqpFunc* func)
{
  return func->type;
}

SLEQP_RETCODE sleqp_func_set_type(SleqpFunc* func,
                                  SLEQP_FUNC_TYPE func_type)
{
  func->type = func_type;

  return SLEQP_OKAY;
}

int sleqp_func_get_num_variables(SleqpFunc* func)
{
  return func->num_variables;
}

int sleqp_func_get_num_constraints(SleqpFunc* func)
{
  return func->num_constraints;
}

SleqpTimer* sleqp_func_get_set_timer(SleqpFunc* func)
{
  return func->set_timer;
}

SleqpTimer* sleqp_func_get_val_timer(SleqpFunc* func)
{
  return func->val_timer;
}

SleqpTimer* sleqp_func_get_grad_timer(SleqpFunc* func)
{
  return func->grad_timer;
}

SleqpTimer* sleqp_func_get_cons_val_timer(SleqpFunc* func)
{
  return func->cons_val_timer;
}

SleqpTimer* sleqp_func_get_cons_jac_timer(SleqpFunc* func)
{
  return func->cons_jac_timer;
}

SleqpTimer* sleqp_func_get_hess_timer(SleqpFunc* func)
{
  return func->hess_timer;
}

SLEQP_RETCODE sleqp_func_hess_prod(SleqpFunc* func,
                                   const double* func_dual,
                                   const SleqpSparseVec* direction,
                                   const SleqpSparseVec* cons_duals,
                                   SleqpSparseVec* product)
{
  assert(func->num_variables == direction->dim);
  assert(func->num_variables == product->dim);
  assert(func->num_constraints == cons_duals->dim);

  assert(sleqp_sparse_vector_is_valid(direction));
  assert(sleqp_sparse_vector_is_finite(direction));

  assert(sleqp_sparse_vector_is_valid(cons_duals));
  assert(sleqp_sparse_vector_is_finite(cons_duals));

  SLEQP_CALL(sleqp_sparse_vector_clear(product));

  SLEQP_CALL(sleqp_timer_start(func->hess_timer));

  SLEQP_CALL(func->callbacks.hess_prod(func,
                                       func_dual,
                                       direction,
                                       cons_duals,
                                       product,
                                       func->data));

  SLEQP_CALL(sleqp_timer_stop(func->hess_timer));

  sleqp_assert_msg(sleqp_sparse_vector_is_valid(product),
                   "Returned invalid Hessian product");

  sleqp_assert_msg(sleqp_sparse_vector_is_finite(product),
                   "Returned Hessian product is not all-finite");

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_func_hess_bilinear(SleqpFunc* func,
                                       const double* func_dual,
                                       const SleqpSparseVec* direction,
                                       const SleqpSparseVec* cons_duals,
                                       double* bilinear_prod)
{
  SLEQP_CALL(sleqp_sparse_vector_clear(func->product));

  SLEQP_CALL(sleqp_func_hess_prod(func,
                                  func_dual,
                                  direction,
                                  cons_duals,
                                  func->product));

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

  SLEQP_CALL(sleqp_timer_free(&func->cons_jac_timer));
  SLEQP_CALL(sleqp_timer_free(&func->cons_val_timer));

  SLEQP_CALL(sleqp_timer_free(&func->grad_timer));
  SLEQP_CALL(sleqp_timer_free(&func->val_timer));
  SLEQP_CALL(sleqp_timer_free(&func->set_timer));

  SLEQP_CALL(sleqp_sparse_vector_free(&func->product));

  SLEQP_CALL(sleqp_hessian_struct_release(&func->hess_struct));

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
