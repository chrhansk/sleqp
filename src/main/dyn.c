#include "dyn.h"

#include "func.h"
#include "mem.h"

typedef struct {
  SleqpDynFuncCallbacks callbacks;
  double accuracy;
  void* func_data;
} DynFuncData;

SLEQP_RETCODE dyn_func_set_value(SleqpFunc* func,
                                 SleqpSparseVec* value,
                                 SLEQP_VALUE_REASON reason,
                                 int* func_grad_nnz,
                                 int* cons_val_nnz,
                                 int* cons_jac_nnz,
                                 void* func_data)
{
  DynFuncData* data = (DynFuncData*) func_data;

  return data->callbacks.set_value(func,
                                   value,
                                   reason,
                                   func_grad_nnz,
                                   cons_val_nnz,
                                   cons_jac_nnz,
                                   data->func_data);
}

SLEQP_RETCODE dyn_func_val(SleqpFunc* func,
                           double* func_val,
                           void* func_data)
{
  DynFuncData* data = (DynFuncData*) func_data;

  return data->callbacks.func_val(func,
                                  data->accuracy,
                                  func_val,
                                  data->func_data);
}

SLEQP_RETCODE dyn_func_grad(SleqpFunc* func,
                            SleqpSparseVec* func_grad,
                            void* func_data)
{
  DynFuncData* data = (DynFuncData*) func_data;

  return data->callbacks.func_grad(func,
                                   func_grad,
                                   data->func_data);
}

SLEQP_RETCODE dyn_func_cons_val(SleqpFunc* func,
                                const SleqpSparseVec* cons_indices,
                                SleqpSparseVec* cons_val,
                                void* func_data)
{
  DynFuncData* data = (DynFuncData*) func_data;

  return data->callbacks.cons_val(func,
                                  data->accuracy,
                                  cons_indices,
                                  cons_val,
                                  data->func_data);
}

SLEQP_RETCODE dyn_func_cons_jac(SleqpFunc* func,
                                const SleqpSparseVec* cons_indices,
                                SleqpSparseMatrix* cons_jac,
                                void* func_data)
{
  DynFuncData* data = (DynFuncData*) func_data;

  return data->callbacks.cons_jac(func,
                                  cons_indices,
                                  cons_jac,
                                  data->func_data);
}

SLEQP_RETCODE dyn_func_hess_product(SleqpFunc* func,
                                    const double* func_dual,
                                    const SleqpSparseVec* direction,
                                    const SleqpSparseVec* cons_duals,
                                    SleqpSparseVec* product,
                                    void* func_data)
{
  DynFuncData* data = (DynFuncData*) func_data;

  return data->callbacks.hess_prod(func,
                                   func_dual,
                                   direction,
                                   cons_duals,
                                   product,
                                   data->func_data);
}

SLEQP_RETCODE dyn_func_free(void* func_data)
{
  if(!func_data)
  {
    return SLEQP_OKAY;
  }

  DynFuncData* data = (DynFuncData*) func_data;

  if(data->callbacks.func_free)
  {
    SLEQP_CALL(data->callbacks.func_free(data->func_data));
  }

  sleqp_free(&data);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_dyn_func_create(SleqpFunc** fstar,
                                    SleqpDynFuncCallbacks* callbacks,
                                    int num_variables,
                                    int num_constraints,
                                    void* func_data)
{
  DynFuncData* data = NULL;

  SLEQP_CALL(sleqp_malloc(&data));

  *data = (DynFuncData) {0};

  data->callbacks = *callbacks;
  data->accuracy = 1.;
  data->func_data = func_data;

  SleqpFuncCallbacks func_callbacks = {
    .set_value = dyn_func_set_value,
    .func_val  = dyn_func_val,
    .func_grad = dyn_func_grad,
    .cons_val  = dyn_func_cons_val,
    .cons_jac  = dyn_func_cons_jac,
    .hess_prod = dyn_func_hess_product,
    .func_free = dyn_func_free
  };

  SLEQP_CALL(sleqp_func_create(fstar,
                               &func_callbacks,
                               num_variables,
                               num_constraints,
                               data));

  SleqpFunc* func = *fstar;

  SLEQP_CALL(sleqp_func_set_type(func, SLEQP_FUNC_TYPE_DYNAMIC));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_dyn_func_set_accuracy(SleqpFunc* func,
                                          double accuracy)
{
  assert(sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_DYNAMIC);

  void* func_data = sleqp_func_get_data(func);

  DynFuncData* data = (DynFuncData*) func_data;

  data->accuracy = accuracy;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_dyn_func_get_accuracy(SleqpFunc* func,
                                          double* accuracy)
{
  assert(sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_DYNAMIC);

  void* func_data = sleqp_func_get_data(func);

  DynFuncData* data = (DynFuncData*) func_data;

  *accuracy = data->accuracy;

  return SLEQP_OKAY;
}
