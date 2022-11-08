#include "func.h"

#include "cmp.h"
#include "fail.h"
#include "log.h"
#include "mem.h"

#include "sparse/mat.h"

struct SleqpFunc
{
  int refcount;

  SleqpFuncCallbacks callbacks;

  SLEQP_FUNC_FLAGS flags;

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

  SleqpVec* product;

  SleqpHessStruct* hess_struct;
};

SLEQP_RETCODE
sleqp_func_create(SleqpFunc** fstar,
                  SleqpFuncCallbacks* callbacks,
                  int num_variables,
                  int num_constraints,
                  void* func_data)
{
  SLEQP_CALL(sleqp_malloc(fstar));

  SleqpFunc* func = *fstar;

  *func          = (SleqpFunc){0};
  func->refcount = 1;

  func->callbacks = *callbacks;

  func->flags = 0;

  func->num_variables   = num_variables;
  func->num_constraints = num_constraints;
  func->data            = func_data;
  func->type            = SLEQP_FUNC_TYPE_REGULAR;

  SLEQP_CALL(sleqp_timer_create(&func->set_timer));
  SLEQP_CALL(sleqp_timer_create(&func->val_timer));
  SLEQP_CALL(sleqp_timer_create(&func->grad_timer));

  SLEQP_CALL(sleqp_timer_create(&func->cons_val_timer));
  SLEQP_CALL(sleqp_timer_create(&func->cons_jac_timer));

  SLEQP_CALL(sleqp_timer_create(&func->hess_timer));

  SLEQP_CALL(sleqp_vec_create_empty(&func->product, num_variables));

  SLEQP_CALL(
    sleqp_hess_struct_create(&func->hess_struct, num_variables, false));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_func_set_value(SleqpFunc* func,
                     SleqpVec* x,
                     SLEQP_VALUE_REASON reason,
                     bool* reject)
{
  assert(sleqp_vec_is_valid(x));
  assert(sleqp_vec_is_finite(x));

  *reject = false;

  SLEQP_CALL(sleqp_timer_start(func->set_timer));

  SLEQP_FUNC_CALL(
    func->callbacks.set_value(func, x, reason, reject, func->data),
    sleqp_func_has_flags(func, SLEQP_FUNC_INTERNAL),
    SLEQP_FUNC_ERROR_SET_VALUE);

  SLEQP_CALL(sleqp_timer_stop(func->set_timer));

  return SLEQP_OKAY;
}

SLEQP_NODISCARD SLEQP_RETCODE
sleqp_func_nonzeros(SleqpFunc* func,
                    int* obj_grad_nnz,
                    int* cons_val_nnz,
                    int* cons_jac_nnz,
                    int* hess_prod_nnz)
{
  *obj_grad_nnz  = SLEQP_NONE;
  *cons_val_nnz  = SLEQP_NONE;
  *cons_jac_nnz  = SLEQP_NONE;
  *hess_prod_nnz = SLEQP_NONE;

  if (func->callbacks.nonzeros)
  {
    SLEQP_FUNC_CALL(func->callbacks.nonzeros(func,
                                             obj_grad_nnz,
                                             cons_val_nnz,
                                             cons_jac_nnz,
                                             hess_prod_nnz,
                                             func->data),
                    sleqp_func_has_flags(func, SLEQP_FUNC_INTERNAL),
                    SLEQP_FUNC_ERROR_NONZEROS);
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_func_obj_val(SleqpFunc* func, double* obj_val)
{
  if (obj_val)
  {
    SLEQP_CALL(sleqp_timer_start(func->val_timer));

    SLEQP_FUNC_CALL(func->callbacks.obj_val(func, obj_val, func->data),
                    sleqp_func_has_flags(func, SLEQP_FUNC_INTERNAL),
                    SLEQP_FUNC_ERROR_OBJ_VAL);

    SLEQP_CALL(sleqp_timer_stop(func->val_timer));

    sleqp_assert_msg(sleqp_is_finite(*obj_val),
                     "Returned infinite function value");
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_func_obj_grad(SleqpFunc* func, SleqpVec* obj_grad)
{
  if (obj_grad)
  {
    SLEQP_CALL(sleqp_vec_clear(obj_grad));

    SLEQP_CALL(sleqp_timer_start(func->grad_timer));

    SLEQP_FUNC_CALL(func->callbacks.obj_grad(func, obj_grad, func->data),
                    sleqp_func_has_flags(func, SLEQP_FUNC_INTERNAL),
                    SLEQP_FUNC_ERROR_OBJ_GRAD);

    SLEQP_CALL(sleqp_timer_stop(func->grad_timer));

    sleqp_assert_msg(sleqp_vec_is_valid(obj_grad),
                     "Returned invalid function gradient");

    sleqp_assert_msg(sleqp_vec_is_finite(obj_grad),
                     "Returned function gradient is not all-finite");
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_func_cons_val(SleqpFunc* func, SleqpVec* cons_val)
{
  const int num_constraints = sleqp_func_num_cons(func);

  if (cons_val)
  {
    SLEQP_CALL(sleqp_vec_clear(cons_val));

    if ((num_constraints != 0))
    {
      SLEQP_CALL(sleqp_timer_start(func->cons_val_timer));

      SLEQP_FUNC_CALL(func->callbacks.cons_val(func, cons_val, func->data),
                      sleqp_func_has_flags(func, SLEQP_FUNC_INTERNAL),
                      SLEQP_FUNC_ERROR_CONS_VAL);

      SLEQP_CALL(sleqp_timer_stop(func->cons_val_timer));
    }

    sleqp_assert_msg(sleqp_vec_is_valid(cons_val),
                     "Returned invalid constraint values");

    sleqp_assert_msg(sleqp_vec_is_finite(cons_val),
                     "Returned constraint values are not all-finite");
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_func_cons_jac(SleqpFunc* func, SleqpMat* cons_jac)
{
  const int num_constraints = sleqp_func_num_cons(func);

  if (cons_jac)
  {
    SLEQP_CALL(sleqp_mat_clear(cons_jac));

    if ((num_constraints != 0))
    {
      SLEQP_CALL(sleqp_timer_start(func->cons_jac_timer));

      SLEQP_FUNC_CALL(func->callbacks.cons_jac(func, cons_jac, func->data),
                      sleqp_func_has_flags(func, SLEQP_FUNC_INTERNAL),
                      SLEQP_FUNC_ERROR_CONS_JAC);

      SLEQP_CALL(sleqp_timer_stop(func->cons_jac_timer));
    }

    sleqp_assert_msg(sleqp_mat_is_valid(cons_jac),
                     "Returned invalid constraint Jacobian");

    sleqp_assert_msg(sleqp_mat_is_finite(cons_jac),
                     "Returned constraint Jacobian is not all-finite");
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_func_eval(SleqpFunc* func,
                double* obj_val,
                SleqpVec* obj_grad,
                SleqpVec* cons_val,
                SleqpMat* cons_jac)
{
  SLEQP_CALL(sleqp_func_obj_val(func, obj_val));

  SLEQP_CALL(sleqp_func_obj_grad(func, obj_grad));

  SLEQP_CALL(sleqp_func_cons_val(func, cons_val));

  SLEQP_CALL(sleqp_func_cons_jac(func, cons_jac));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_func_set_callbacks(SleqpFunc* func, SleqpFuncCallbacks* callbacks)
{
  func->callbacks = *callbacks;

  return SLEQP_OKAY;
}

SleqpHessStruct*
sleqp_func_hess_struct(SleqpFunc* func)
{
  return func->hess_struct;
}

SLEQP_RETCODE
sleqp_func_flags_set(SleqpFunc* func, SLEQP_FUNC_FLAGS flags, bool value)
{
  if (value)
  {
    return sleqp_func_flags_add(func, flags);
  }
  else
  {
    return sleqp_func_flags_remove(func, flags);
  }
}

SLEQP_RETCODE
sleqp_func_flags_add(SleqpFunc* func, SLEQP_FUNC_FLAGS flags)
{
  func->flags |= flags;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_func_flags_remove(SleqpFunc* func, SLEQP_FUNC_FLAGS flags)
{
  func->flags &= ~(flags);

  return SLEQP_OKAY;
}

bool
sleqp_func_has_flags(const SleqpFunc* func, SLEQP_FUNC_FLAGS flags)
{
  return !!(func->flags & flags);
}

bool
sleqp_func_flags_copy(const SleqpFunc* source,
                      SleqpFunc* target,
                      SLEQP_FUNC_FLAGS flags)
{
  SLEQP_CALL(
    sleqp_func_flags_set(target, flags, sleqp_func_has_flags(source, flags)));
  return SLEQP_OKAY;
}

SLEQP_FUNC_TYPE
sleqp_func_get_type(const SleqpFunc* func)
{
  return func->type;
}

SLEQP_RETCODE
sleqp_func_set_type(SleqpFunc* func, SLEQP_FUNC_TYPE func_type)
{
  func->type = func_type;

  return SLEQP_OKAY;
}

int
sleqp_func_num_vars(const SleqpFunc* func)
{
  return func->num_variables;
}

int
sleqp_func_num_cons(const SleqpFunc* func)
{
  return func->num_constraints;
}

SleqpTimer*
sleqp_func_get_set_timer(SleqpFunc* func)
{
  return func->set_timer;
}

SleqpTimer*
sleqp_func_get_val_timer(SleqpFunc* func)
{
  return func->val_timer;
}

SleqpTimer*
sleqp_func_get_grad_timer(SleqpFunc* func)
{
  return func->grad_timer;
}

SleqpTimer*
sleqp_func_get_cons_val_timer(SleqpFunc* func)
{
  return func->cons_val_timer;
}

SleqpTimer*
sleqp_func_get_cons_jac_timer(SleqpFunc* func)
{
  return func->cons_jac_timer;
}

SleqpTimer*
sleqp_func_get_hess_timer(SleqpFunc* func)
{
  return func->hess_timer;
}

SLEQP_RETCODE
sleqp_func_hess_prod(SleqpFunc* func,
                     const double* obj_dual,
                     const SleqpVec* direction,
                     const SleqpVec* cons_duals,
                     SleqpVec* product)
{
  assert(func->num_variables == direction->dim);
  assert(func->num_variables == product->dim);
  assert(func->num_constraints == cons_duals->dim);

  assert(sleqp_vec_is_valid(direction));
  assert(sleqp_vec_is_finite(direction));

  assert(sleqp_vec_is_valid(cons_duals));
  assert(sleqp_vec_is_finite(cons_duals));

  SLEQP_CALL(sleqp_vec_clear(product));

  SLEQP_CALL(sleqp_timer_start(func->hess_timer));

  SLEQP_FUNC_CALL(
    func->callbacks
      .hess_prod(func, obj_dual, direction, cons_duals, product, func->data),
    sleqp_func_has_flags(func, SLEQP_FUNC_INTERNAL | SLEQP_FUNC_HESS_INTERNAL),
    SLEQP_FUNC_ERROR_HESS_PROD);

  SLEQP_CALL(sleqp_timer_stop(func->hess_timer));

  sleqp_assert_msg(sleqp_vec_is_valid(product),
                   "Returned invalid Hessian product");

  sleqp_assert_msg(sleqp_vec_is_finite(product),
                   "Returned Hessian product is not all-finite");

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_func_hess_bilinear(SleqpFunc* func,
                         const double* obj_dual,
                         const SleqpVec* direction,
                         const SleqpVec* cons_duals,
                         double* bilinear_prod)
{
  SLEQP_CALL(sleqp_vec_clear(func->product));

  SLEQP_CALL(
    sleqp_func_hess_prod(func, obj_dual, direction, cons_duals, func->product));

  SLEQP_CALL(sleqp_vec_dot(direction, func->product, bilinear_prod));

  return SLEQP_OKAY;
}

void*
sleqp_func_get_data(SleqpFunc* func)
{
  return func->data;
}

static SLEQP_RETCODE
func_free(SleqpFunc** fstar)
{
  SleqpFunc* func = *fstar;

  if (!func)
  {
    return SLEQP_OKAY;
  }

  if (func->callbacks.func_free)
  {
    SLEQP_CALL(func->callbacks.func_free(func->data));
  }

  SLEQP_CALL(sleqp_timer_free(&func->hess_timer));

  SLEQP_CALL(sleqp_timer_free(&func->cons_jac_timer));
  SLEQP_CALL(sleqp_timer_free(&func->cons_val_timer));

  SLEQP_CALL(sleqp_timer_free(&func->grad_timer));
  SLEQP_CALL(sleqp_timer_free(&func->val_timer));
  SLEQP_CALL(sleqp_timer_free(&func->set_timer));

  SLEQP_CALL(sleqp_vec_free(&func->product));

  SLEQP_CALL(sleqp_hess_struct_release(&func->hess_struct));

  sleqp_free(fstar);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_func_capture(SleqpFunc* func)
{
  ++func->refcount;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_func_release(SleqpFunc** fstar)
{
  SleqpFunc* func = *fstar;

  if (!func)
  {
    return SLEQP_OKAY;
  }

  if (--func->refcount == 0)
  {
    SLEQP_CALL(func_free(fstar));
  }

  *fstar = NULL;

  return SLEQP_OKAY;
}
