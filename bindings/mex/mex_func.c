#include "mex_func.h"

#include <assert.h>

#include "mex_error.h"
#include "mex_fields.h"
#include "mex_func_common.h"
#include "mex_hess.h"
#include "sleqp/pub_types.h"

/*
static void
print_array(mxArray* array)
{
  mexCallMATLAB(0, NULL, 1, &array, MATLAB_FUNC_DISP);
}
*/

// TODO: Better handling of matlab exceptions
// TODO: Better assert

typedef struct
{
  SleqpSettings* settings;

  // Callbacks
  struct
  {
    mxArray* obj_val;
    mxArray* obj_grad;
    mxArray* cons_val;
    mxArray* cons_jac;
  } callbacks;

  mxArray* primal;

  MexHess hess;

} FuncData;

static SLEQP_RETCODE
mex_func_set(SleqpFunc* func,
             SleqpVec* primal,
             SLEQP_VALUE_REASON reason,
             bool* reject,
             void* data)
{
  FuncData* func_data = (FuncData*)data;

  SLEQP_CALL(sleqp_vec_to_raw(primal, mxGetPr(func_data->primal)));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
mex_func_obj_val(SleqpFunc* func, double* obj_val, void* data)
{
  FuncData* func_data = (FuncData*)data;

  mxArray* rhs[] = {func_data->callbacks.obj_val, func_data->primal};

  MEX_EVAL_INTO_REAL(rhs, obj_val);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
mex_func_obj_grad(SleqpFunc* func, SleqpVec* obj_grad, void* data)
{
  FuncData* func_data = (FuncData*)data;

  mxArray* rhs[] = {func_data->callbacks.obj_grad, func_data->primal};

  MEX_EVAL_INTO_VEC(rhs, func_data->settings, obj_grad);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
mex_func_cons_val(SleqpFunc* func, SleqpVec* cons_val, void* data)
{
  FuncData* func_data = (FuncData*)data;

  mxArray* rhs[] = {func_data->callbacks.cons_val, func_data->primal};

  MEX_EVAL_INTO_VEC(rhs, func_data->settings, cons_val);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
mex_func_cons_jac(SleqpFunc* func, SleqpMat* cons_jac, void* data)
{
  FuncData* func_data = (FuncData*)data;

  mxArray* rhs[] = {func_data->callbacks.cons_jac, func_data->primal};

  MEX_EVAL_INTO_SPARSE_MATRIX(rhs, func_data->settings, cons_jac);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
mex_func_hess_prod(SleqpFunc* func,
                   const SleqpVec* direction,
                   const SleqpVec* cons_duals,
                   SleqpVec* result,
                   void* data)
{
  FuncData* func_data = (FuncData*)data;

  SLEQP_CALL(mex_hess_prod(&(func_data->hess),
                           func_data->primal,
                           direction,
                           cons_duals,
                           NULL,
                           0,
                           result));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
create_func_data(FuncData** star,
                 const mxArray* mex_callbacks,
                 SleqpSettings* settings,
                 int num_vars,
                 int num_cons,
                 bool with_hess)
{
  SLEQP_CALL(sleqp_malloc(star));

  FuncData* func_data = *star;

  *func_data = (FuncData){0};

  SLEQP_CALL(sleqp_settings_capture(settings));
  func_data->settings = settings;

  assert(mxIsStruct(mex_callbacks));

  SLEQP_CALL(mex_callback_from_struct(mex_callbacks,
                                      MEX_INPUT_OBJ_VAL,
                                      &func_data->callbacks.obj_val));

  SLEQP_CALL(mex_callback_from_struct(mex_callbacks,
                                      MEX_INPUT_OBJ_GRAD,
                                      &func_data->callbacks.obj_grad));

  if (num_cons > 0)
  {
    SLEQP_CALL(mex_callback_from_struct(mex_callbacks,
                                        MEX_INPUT_CONS_VAL,
                                        &func_data->callbacks.cons_val));

    SLEQP_CALL(mex_callback_from_struct(mex_callbacks,
                                        MEX_INPUT_CONS_JAC,
                                        &func_data->callbacks.cons_jac));
  }

  if (with_hess)
  {
    SLEQP_CALL(mex_hess_init(&(func_data->hess),
                             settings,
                             mex_callbacks,
                             num_vars,
                             num_cons));
  }

  func_data->primal = mxCreateDoubleMatrix(1, num_vars, mxREAL);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
mex_func_free(void* data)
{
  FuncData* func_data = (FuncData*)data;

  SLEQP_CALL(mex_hess_free(&func_data->hess));

  mxDestroyArray(func_data->primal);

  SLEQP_CALL(sleqp_settings_release(&func_data->settings));

  sleqp_free(&func_data);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
mex_func_create(SleqpFunc** star,
                const mxArray* mex_callbacks,
                SleqpSettings* settings,
                int num_vars,
                int num_cons,
                bool with_hess)
{
  SleqpFuncCallbacks callbacks = {.set_value = mex_func_set,
                                  .obj_val   = mex_func_obj_val,
                                  .obj_grad  = mex_func_obj_grad,
                                  .cons_val  = mex_func_cons_val,
                                  .cons_jac  = mex_func_cons_jac,
                                  .hess_prod = mex_func_hess_prod,
                                  .func_free = mex_func_free};

  FuncData* func_data;

  SLEQP_CALL(create_func_data(&func_data,
                              mex_callbacks,
                              settings,
                              num_vars,
                              num_cons,
                              with_hess));

  SLEQP_CALL(
    sleqp_func_create(star, &callbacks, num_vars, num_cons, (void*)func_data));

  return SLEQP_OKAY;
}
