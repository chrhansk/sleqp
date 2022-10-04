#include "dyn.h"

#include "func.h"
#include "mem.h"
#include "pub_func.h"
#include "sparse/vec.h"

typedef enum
{
  HAS_OBJ_WEIGHT   = (1 << 0),
  HAS_CONS_WEIGHTS = (1 << 1),
  HAS_ERROR_BOUND  = (1 << 2),
  HAS_VALUES       = (1 << 3)
} FLAGS;

typedef struct
{
  SleqpDynFuncCallbacks callbacks;
  void* func_data;
  FLAGS flags;
  double obj_val;
  SleqpVec* cons_val;
  double error_bound;
  double error;

} DynFuncData;

#define ERROR_SET_ERROR_BOUND "Error '%s' setting error bound"
#define ERROR_SET_OBJ_WEIGHT "Error '%s' setting objective weight"
#define ERROR_SET_CONS_WEIGHTS "Error '%s' setting constraint weights"

static SLEQP_RETCODE
dyn_func_set_value(SleqpFunc* func,
                   SleqpVec* value,
                   SLEQP_VALUE_REASON reason,
                   bool* reject,
                   void* func_data)
{
  DynFuncData* data = (DynFuncData*)func_data;

  data->flags &= ~(HAS_VALUES);
  data->error = SLEQP_NONE;

  SLEQP_FUNC_CALL(
    data->callbacks.set_value(func, value, reason, reject, data->func_data),
    sleqp_func_has_flags(func, SLEQP_FUNC_INTERNAL),
    SLEQP_FUNC_ERROR_SET_VALUE);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
dyn_func_nonzeros(SleqpFunc* func,
                  int* obj_grad_nnz,
                  int* cons_val_nnz,
                  int* cons_jac_nnz,
                  int* hess_prod_nnz,
                  void* func_data)
{
  *obj_grad_nnz  = SLEQP_NONE;
  *cons_val_nnz  = SLEQP_NONE;
  *cons_jac_nnz  = SLEQP_NONE;
  *hess_prod_nnz = SLEQP_NONE;

  DynFuncData* data = (DynFuncData*)func_data;

  if (data->callbacks.nonzeros)
  {
    SLEQP_FUNC_CALL(data->callbacks.nonzeros(func,
                                             obj_grad_nnz,
                                             cons_val_nnz,
                                             cons_jac_nnz,
                                             hess_prod_nnz,
                                             data->func_data),
                    sleqp_func_has_flags(func, SLEQP_FUNC_INTERNAL),
                    SLEQP_FUNC_ERROR_NONZEROS);

    if (*cons_val_nnz != SLEQP_NONE)
    {
      SLEQP_CALL(sleqp_vec_reserve(data->cons_val, *cons_val_nnz));
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
dyn_func_eval(SleqpFunc* func, DynFuncData* data)
{
#if SLEQP_DEBUG
  {
    const int num_cons = sleqp_func_num_cons(func);

    if (num_cons != 0)
    {
      assert(data->flags & HAS_CONS_WEIGHTS);
    }

    assert(data->flags & HAS_OBJ_WEIGHT);
    assert(data->flags & HAS_ERROR_BOUND);
  }
#endif

  if (data->flags & HAS_VALUES)
  {
    assert(data->error != SLEQP_NONE);
    return SLEQP_OKAY;
  }

  data->error = SLEQP_NONE;

  SLEQP_FUNC_CALL(data->callbacks.eval(func,
                                       &data->obj_val,
                                       data->cons_val,
                                       &data->error,
                                       data->func_data),
                  sleqp_func_has_flags(func, SLEQP_FUNC_INTERNAL),
                  SLEQP_FUNC_ERROR_OBJ_GRAD);

  if (data->error == SLEQP_NONE)
  {
    sleqp_raise(SLEQP_FUNC_EVAL_ERROR, "Dynamic function did not report error");
  }
  else if (data->error < 0)
  {
    sleqp_raise(SLEQP_FUNC_EVAL_ERROR,
                "Dynamic function reported negative error of %e",
                data->error);
  }
  else if (data->error > data->error_bound)
  {
    sleqp_raise(SLEQP_FUNC_EVAL_ERROR,
                "Dynamic function error of %e exceeds error bound of %e",
                data->error,
                data->error_bound);
  }

  data->flags |= HAS_VALUES;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
dyn_func_obj_val(SleqpFunc* func, double* obj_val, void* func_data)
{
  DynFuncData* data = (DynFuncData*)func_data;

  SLEQP_CALL(dyn_func_eval(func, data));

  *obj_val = data->obj_val;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
dyn_func_obj_grad(SleqpFunc* func, SleqpVec* obj_grad, void* func_data)
{
  DynFuncData* data = (DynFuncData*)func_data;

  SLEQP_FUNC_CALL(data->callbacks.obj_grad(func, obj_grad, data->func_data),
                  sleqp_func_has_flags(func, SLEQP_FUNC_INTERNAL),
                  SLEQP_FUNC_ERROR_OBJ_GRAD);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
dyn_func_cons_val(SleqpFunc* func, SleqpVec* cons_val, void* func_data)
{
  DynFuncData* data = (DynFuncData*)func_data;

  SLEQP_CALL(dyn_func_eval(func, data));

  SLEQP_CALL(sleqp_vec_copy(data->cons_val, cons_val));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
dyn_func_cons_jac(SleqpFunc* func, SleqpSparseMatrix* cons_jac, void* func_data)
{
  DynFuncData* data = (DynFuncData*)func_data;

  SLEQP_FUNC_CALL(data->callbacks.cons_jac(func, cons_jac, data->func_data),
                  sleqp_func_has_flags(func, SLEQP_FUNC_INTERNAL),
                  SLEQP_FUNC_ERROR_CONS_JAC);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
dyn_func_hess_product(SleqpFunc* func,
                      const double* obj_dual,
                      const SleqpVec* direction,
                      const SleqpVec* cons_duals,
                      SleqpVec* product,
                      void* func_data)
{
  DynFuncData* data = (DynFuncData*)func_data;

  SLEQP_FUNC_CALL(
    data->callbacks.hess_prod(func,
                              obj_dual,
                              direction,
                              cons_duals,
                              product,
                              data->func_data),
    sleqp_func_has_flags(func, SLEQP_FUNC_INTERNAL | SLEQP_FUNC_HESS_INTERNAL),
    SLEQP_FUNC_ERROR_HESS_PROD);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
dyn_func_free(void* func_data)
{
  if (!func_data)
  {
    return SLEQP_OKAY;
  }

  DynFuncData* data = (DynFuncData*)func_data;

  if (data->callbacks.func_free)
  {
    SLEQP_CALL(data->callbacks.func_free(data->func_data));
  }

  SLEQP_CALL(sleqp_vec_free(&data->cons_val));

  sleqp_free(&data);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_dyn_func_create(SleqpFunc** fstar,
                      SleqpDynFuncCallbacks* callbacks,
                      int num_variables,
                      int num_constraints,
                      void* func_data)
{
  DynFuncData* data = NULL;

  SLEQP_CALL(sleqp_malloc(&data));

  *data = (DynFuncData){0};

  data->callbacks = *callbacks;
  data->func_data = func_data;
  data->flags     = 0;

  SleqpFuncCallbacks func_callbacks = {.set_value = dyn_func_set_value,
                                       .nonzeros  = dyn_func_nonzeros,
                                       .obj_val   = dyn_func_obj_val,
                                       .obj_grad  = dyn_func_obj_grad,
                                       .cons_val  = dyn_func_cons_val,
                                       .cons_jac  = dyn_func_cons_jac,
                                       .hess_prod = dyn_func_hess_product,
                                       .func_free = dyn_func_free};

  SLEQP_CALL(sleqp_func_create(fstar,
                               &func_callbacks,
                               num_variables,
                               num_constraints,
                               data));

  SLEQP_CALL(sleqp_vec_create_empty(&data->cons_val, num_constraints));

  SleqpFunc* func = *fstar;

  SLEQP_CALL(sleqp_func_set_type(func, SLEQP_FUNC_TYPE_DYNAMIC));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_dyn_func_eval(SleqpFunc* func,
                    double* obj_val,
                    SleqpVec* cons_val,
                    double* error)
{
  DynFuncData* data = (DynFuncData*)sleqp_func_get_data(func);

  SLEQP_CALL(dyn_func_eval(func, data));

  *obj_val = data->obj_val;
  *error   = data->error;
  SLEQP_CALL(sleqp_vec_copy(data->cons_val, cons_val));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_dyn_func_set_error_bound(SleqpFunc* func, double error_bound)
{
  assert(sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_DYNAMIC);

  DynFuncData* data = (DynFuncData*)sleqp_func_get_data(func);

  data->flags |= HAS_ERROR_BOUND;
  data->flags &= ~(HAS_VALUES);
  data->error_bound = error_bound;

  SLEQP_FUNC_CALL(
    data->callbacks.set_error_bound(func, error_bound, data->func_data),
    sleqp_func_has_flags(func, SLEQP_FUNC_INTERNAL),
    ERROR_SET_ERROR_BOUND);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_dyn_func_error(SleqpFunc* func, double* error)
{
  assert(sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_DYNAMIC);

  DynFuncData* data = (DynFuncData*)sleqp_func_get_data(func);

  *error = data->error;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_dyn_func_error_bound(SleqpFunc* func, double* error_bound)
{
  assert(sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_DYNAMIC);

  DynFuncData* data = (DynFuncData*)sleqp_func_get_data(func);

  *error_bound = data->error_bound;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_dyn_func_error_estimate(SleqpFunc* func, double* error_estimate)
{
  assert(sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_DYNAMIC);

  DynFuncData* data = (DynFuncData*)sleqp_func_get_data(func);

  if (data->error != SLEQP_NONE)
  {
    *error_estimate = data->error;
  }
  else
  {
    *error_estimate = data->error_bound;
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_dyn_func_set_obj_weight(SleqpFunc* func, double obj_weight)
{
  assert(sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_DYNAMIC);

  DynFuncData* data = (DynFuncData*)sleqp_func_get_data(func);

  data->flags |= HAS_OBJ_WEIGHT;
  data->flags &= ~(1 << HAS_VALUES);

  SLEQP_FUNC_CALL(
    data->callbacks.set_obj_weight(func, obj_weight, data->func_data),
    sleqp_func_has_flags(func, SLEQP_FUNC_INTERNAL),
    ERROR_SET_OBJ_WEIGHT);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_dyn_func_set_cons_weights(SleqpFunc* func, const double* cons_weights)
{
  assert(sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_DYNAMIC);

  DynFuncData* data = (DynFuncData*)sleqp_func_get_data(func);

  data->flags |= HAS_CONS_WEIGHTS;
  data->flags &= ~(1 << HAS_VALUES);

  const int num_cons = sleqp_func_num_cons(func);

  if (num_cons != 0)
  {
    SLEQP_FUNC_CALL(
      data->callbacks.set_cons_weights(func, cons_weights, data->func_data),
      sleqp_func_has_flags(func, SLEQP_FUNC_INTERNAL),
      ERROR_SET_CONS_WEIGHTS);
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_dyn_func_set_callbacks(SleqpFunc* func, SleqpDynFuncCallbacks* callbacks)
{
  assert(sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_DYNAMIC);
  DynFuncData* data = (DynFuncData*)sleqp_func_get_data(func);
  assert(data);

  data->callbacks = *callbacks;

  return SLEQP_OKAY;
}
