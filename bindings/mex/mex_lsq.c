#include "mex_lsq.h"

#include "mex_fields.h"
#include "mex_func_common.h"

typedef struct
{
  SleqpSettings* settings;

  // Callbacks
  struct
  {
    mxArray* lsq_residuals;
    mxArray* lsq_jac_forward;
    mxArray* lsq_jac_adjoint;
    mxArray* cons_val;
    mxArray* cons_jac;
  } callbacks;

  mxArray* primal;

  mxArray* forward_dir;
  mxArray* adjoint_dir;

} LSQFuncData;

SLEQP_RETCODE
get_num_residuals(const mxArray* mex_x0,
                  const mxArray* mex_res_callback,
                  int* num_residuals)
{
  mxArray* lhs[] = {NULL};

  // Need to strip const'ness away
  // due to signature of mexCallMATLABWithTrap
  mxArray* rhs[] = {(mxArray*)mex_res_callback, (mxArray*)mex_x0};

  MEX_EVAL(lhs, rhs);

  *num_residuals = mxGetNumberOfElements(*lhs);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
create_lsq_func_data(LSQFuncData** star,
                     const mxArray* mex_x0,
                     const mxArray* mex_callbacks,
                     SleqpSettings* settings,
                     int num_variables,
                     int num_constraints,
                     int* num_residuals)
{
  SLEQP_CALL(sleqp_malloc(star));

  LSQFuncData* func_data = *star;

  *func_data = (LSQFuncData){0};

  SLEQP_CALL(sleqp_settings_capture(settings));
  func_data->settings = settings;

  SLEQP_CALL(mex_callback_from_struct(mex_callbacks,
                                      MEX_INPUT_LSQ_RES,
                                      &func_data->callbacks.lsq_residuals));

  SLEQP_CALL(mex_callback_from_struct(mex_callbacks,
                                      MEX_INPUT_LSQ_JAC_FWD,
                                      &func_data->callbacks.lsq_jac_forward));

  SLEQP_CALL(mex_callback_from_struct(mex_callbacks,
                                      MEX_INPUT_LSQ_JAC_ADJ,
                                      &func_data->callbacks.lsq_jac_adjoint));

  if (num_constraints > 0)
  {
    SLEQP_CALL(mex_callback_from_struct(mex_callbacks,
                                        MEX_INPUT_CONS_VAL,
                                        &func_data->callbacks.cons_val));

    SLEQP_CALL(mex_callback_from_struct(mex_callbacks,
                                        MEX_INPUT_CONS_JAC,
                                        &func_data->callbacks.cons_jac));
  }

  SLEQP_CALL(get_num_residuals(mex_x0,
                               func_data->callbacks.lsq_residuals,
                               num_residuals));

  func_data->primal = mxCreateDoubleMatrix(1, num_variables, mxREAL);

  func_data->adjoint_dir = mxCreateDoubleMatrix(1, *num_residuals, mxREAL);
  func_data->forward_dir = mxCreateDoubleMatrix(1, num_variables, mxREAL);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
mex_lsq_func_set(SleqpFunc* func,
                 SleqpVec* x,
                 SLEQP_VALUE_REASON reason,
                 bool* reject,
                 void* data)
{
  LSQFuncData* func_data = (LSQFuncData*)data;

  SLEQP_CALL(sleqp_vec_to_raw(x, mxGetPr(func_data->primal)));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
mex_lsq_func_cons_val(SleqpFunc* func, SleqpVec* cons_val, void* data)
{
  LSQFuncData* func_data = (LSQFuncData*)data;

  mxArray* rhs[] = {func_data->callbacks.cons_val, func_data->primal};

  MEX_EVAL_INTO_VEC(rhs, func_data->settings, cons_val);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
mex_lsq_func_cons_jac(SleqpFunc* func, SleqpMat* cons_jac, void* data)
{
  LSQFuncData* func_data = (LSQFuncData*)data;

  mxArray* rhs[] = {func_data->callbacks.cons_jac, func_data->primal};

  MEX_EVAL_INTO_SPARSE_MATRIX(rhs, func_data->settings, cons_jac);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
mex_lsq_residuals(SleqpFunc* func, SleqpVec* residual, void* data)
{
  LSQFuncData* func_data = (LSQFuncData*)data;

  mxArray* rhs[] = {func_data->callbacks.lsq_residuals, func_data->primal};

  MEX_EVAL_INTO_VEC(rhs, func_data->settings, residual);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
mex_lsq_jac_forward(SleqpFunc* func,
                    const SleqpVec* forward_direction,
                    SleqpVec* product,
                    void* data)
{
  LSQFuncData* func_data = (LSQFuncData*)data;

  SLEQP_CALL(
    sleqp_vec_to_raw(forward_direction, mxGetPr(func_data->forward_dir)));

  mxArray* rhs[] = {func_data->callbacks.lsq_jac_forward,
                    func_data->primal,
                    func_data->forward_dir};

  MEX_EVAL_INTO_VEC(rhs, func_data->settings, product);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
mex_lsq_jac_adjoint(SleqpFunc* func,
                    const SleqpVec* adjoint_direction,
                    SleqpVec* product,
                    void* data)
{
  LSQFuncData* func_data = (LSQFuncData*)data;

  SLEQP_CALL(
    sleqp_vec_to_raw(adjoint_direction, mxGetPr(func_data->adjoint_dir)));

  mxArray* rhs[] = {func_data->callbacks.lsq_jac_adjoint,
                    func_data->primal,
                    func_data->adjoint_dir};

  MEX_EVAL_INTO_VEC(rhs, func_data->settings, product);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
mex_lsq_func_free(void* data)
{
  LSQFuncData* func_data = (LSQFuncData*)data;

  mxDestroyArray(func_data->adjoint_dir);
  mxDestroyArray(func_data->forward_dir);
  mxDestroyArray(func_data->primal);

  SLEQP_CALL(sleqp_settings_release(&func_data->settings));

  sleqp_free(&func_data);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
mex_lsq_func_create(SleqpFunc** star,
                    const mxArray* mex_x0,
                    const mxArray* mex_callbacks,
                    SleqpSettings* settings,
                    int num_variables,
                    int num_constraints)
{
  SleqpLSQCallbacks callbacks = {.set_value       = mex_lsq_func_set,
                                 .lsq_residuals   = mex_lsq_residuals,
                                 .lsq_jac_forward = mex_lsq_jac_forward,
                                 .lsq_jac_adjoint = mex_lsq_jac_adjoint,
                                 .cons_val        = mex_lsq_func_cons_val,
                                 .cons_jac        = mex_lsq_func_cons_jac,
                                 .func_free       = mex_lsq_func_free};

  int num_residuals;

  LSQFuncData* func_data;

  SLEQP_CALL(create_lsq_func_data(&func_data,
                                  mex_x0,
                                  mex_callbacks,
                                  settings,
                                  num_variables,
                                  num_constraints,
                                  &num_residuals));

  SLEQP_CALL(sleqp_lsq_func_create(star,
                                   &callbacks,
                                   num_variables,
                                   num_constraints,
                                   num_residuals,
                                   0.,
                                   settings,
                                   func_data));

  return SLEQP_OKAY;
}
