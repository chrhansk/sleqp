#include "dyn.h"

#include "func.h"
#include "mem.h"

typedef struct
{
  SleqpDynFuncCallbacks callbacks;
  double accuracy;
  void* func_data;
} DynFuncData;

static SLEQP_RETCODE
dyn_func_set_value(SleqpFunc* func,
                   SleqpVec* value,
                   SLEQP_VALUE_REASON reason,
                   bool* reject,
                   void* func_data)
{
  DynFuncData* data = (DynFuncData*)func_data;

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
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
dyn_func_obj_val(SleqpFunc* func, double* obj_val, void* func_data)
{
  DynFuncData* data = (DynFuncData*)func_data;

  SLEQP_FUNC_CALL(
    data->callbacks.obj_val(func, data->accuracy, obj_val, data->func_data),
    sleqp_func_has_flags(func, SLEQP_FUNC_INTERNAL),
    SLEQP_FUNC_ERROR_OBJ_VAL);

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

  SLEQP_FUNC_CALL(
    data->callbacks.cons_val(func, data->accuracy, cons_val, data->func_data),
    sleqp_func_has_flags(func, SLEQP_FUNC_INTERNAL),
    SLEQP_FUNC_ERROR_CONS_VAL);

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
  data->accuracy  = 1.;
  data->func_data = func_data;

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

  SleqpFunc* func = *fstar;

  SLEQP_CALL(sleqp_func_set_type(func, SLEQP_FUNC_TYPE_DYNAMIC));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_dyn_func_set_accuracy(SleqpFunc* func, double accuracy)
{
  assert(sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_DYNAMIC);

  void* func_data = sleqp_func_get_data(func);

  DynFuncData* data = (DynFuncData*)func_data;

  data->accuracy = accuracy;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_dyn_func_get_accuracy(SleqpFunc* func, double* accuracy)
{
  assert(sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_DYNAMIC);

  void* func_data = sleqp_func_get_data(func);

  DynFuncData* data = (DynFuncData*)func_data;

  *accuracy = data->accuracy;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_dyn_func_obj_val(SleqpFunc* func, double accuracy, double* obj_val)
{
  assert(sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_DYNAMIC);

  void* func_data = sleqp_func_get_data(func);

  DynFuncData* data = (DynFuncData*)func_data;

  SLEQP_FUNC_CALL(
    data->callbacks.obj_val(func, accuracy, obj_val, data->func_data),
    sleqp_func_has_flags(func, SLEQP_FUNC_INTERNAL),
    SLEQP_FUNC_ERROR_OBJ_VAL);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_dyn_func_cons_val(SleqpFunc* func, double accuracy, SleqpVec* cons_val)
{
  assert(sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_DYNAMIC);

  void* func_data = sleqp_func_get_data(func);

  DynFuncData* data = (DynFuncData*)func_data;

  SLEQP_FUNC_CALL(
    data->callbacks.cons_val(func, accuracy, cons_val, data->func_data),
    sleqp_func_has_flags(func, SLEQP_FUNC_INTERNAL),
    SLEQP_FUNC_ERROR_CONS_VAL);

  return SLEQP_OKAY;
}
