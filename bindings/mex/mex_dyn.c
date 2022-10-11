#include "mex_lsq.h"

#include "mex_error.h"
#include "mex_fields.h"
#include "mex_func_common.h"
#include "mex_hess.h"

typedef struct
{
  SleqpParams* params;

  // Callbacks
  struct
  {
    mxArray* eval;

    mxArray* obj_grad;
    mxArray* cons_jac;
  } callbacks;

  mxArray* primal;

  mxArray* error_bound;
  mxArray* obj_weight;
  mxArray* cons_weights;
  mxArray* error_info;

  MexHess hess;

} DynFuncData;
static const char* field_names[]
  = {MEX_INPUT_DYN_ERROR_BOUND, MEX_INPUT_DYN_OBJ_WEIGHT};

static const char* field_names_cons[] = {MEX_INPUT_DYN_ERROR_BOUND,
                                         MEX_INPUT_DYN_OBJ_WEIGHT,
                                         MEX_INPUT_DYN_CONS_WEIGHTS};

static SLEQP_RETCODE
create_dyn_func_data(DynFuncData** star,
                     const mxArray* mex_callbacks,
                     SleqpParams* params,
                     int num_vars,
                     int num_cons,
                     bool with_hess)
{
  SLEQP_CALL(sleqp_malloc(star));

  DynFuncData* func_data = *star;

  *func_data = (DynFuncData){0};

  SLEQP_CALL(sleqp_params_capture(params));
  func_data->params = params;

  MEX_EXPECT_STRUCT(mex_callbacks);

  SLEQP_CALL(mex_callback_from_struct(mex_callbacks,
                                      MEX_INPUT_DYN_EVAL,
                                      &func_data->callbacks.eval));

  SLEQP_CALL(mex_callback_from_struct(mex_callbacks,
                                      MEX_INPUT_OBJ_GRAD,
                                      &func_data->callbacks.obj_grad));

  if (num_cons > 0)
  {
    SLEQP_CALL(mex_callback_from_struct(mex_callbacks,
                                        MEX_INPUT_CONS_JAC,
                                        &func_data->callbacks.cons_jac));

    func_data->cons_weights = mxCreateDoubleMatrix(1, num_cons, mxREAL);
  }

  if (with_hess)
  {
    SLEQP_CALL(mex_hess_init(&(func_data->hess),
                             params,
                             mex_callbacks,
                             num_vars,
                             num_cons));
  }

  func_data->primal = mxCreateDoubleMatrix(1, num_vars, mxREAL);

  func_data->error_bound = mxCreateDoubleScalar(0.);
  func_data->obj_weight  = mxCreateDoubleScalar(0.);

  if (num_cons > 0)
  {
    func_data->error_info
      = mxCreateStructMatrix(1,
                             1,
                             MEX_ARRAY_LEN(field_names_cons),
                             field_names_cons);

    mxSetFieldByNumber(func_data->error_info, 0, 2, func_data->cons_weights);
  }
  else
  {
    func_data->error_info
      = mxCreateStructMatrix(1, 1, MEX_ARRAY_LEN(field_names), field_names);
  }

  mxSetFieldByNumber(func_data->error_info, 0, 0, func_data->error_bound);
  mxSetFieldByNumber(func_data->error_info, 0, 1, func_data->obj_weight);

  if (num_cons > 0)
  {
    mxSetFieldByNumber(func_data->error_info, 0, 2, func_data->cons_weights);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
dyn_func_set(SleqpFunc* func,
             SleqpVec* primal,
             SLEQP_VALUE_REASON reason,
             bool* reject,
             void* data)
{
  DynFuncData* func_data = (DynFuncData*)data;

  SLEQP_CALL(sleqp_vec_to_raw(primal, mxGetPr(func_data->primal)));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
dyn_func_eval(SleqpFunc* func,
              double* obj_val,
              SleqpVec* cons_val,
              double* error,
              void* data)
{
  DynFuncData* func_data = (DynFuncData*)data;

  mxArray* lhs[] = {NULL, NULL, NULL};

  mxArray* rhs[]
    = {func_data->callbacks.eval, func_data->primal, func_data->error_info};

  MEX_EVAL(lhs, rhs);

  {
    mxArray* obj_array = lhs[0];

    SLEQP_CALL(mex_array_to_real(obj_array, obj_val));
  }

  {
    mxArray* cons_array = lhs[1];

    SLEQP_CALL(mex_array_to_vec(cons_array, func_data->params, cons_val));
  }

  {
    mxArray* error_array = lhs[2];

    SLEQP_CALL(mex_array_to_real(error_array, error));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
dyn_func_obj_grad(SleqpFunc* func, SleqpVec* obj_grad, void* data)
{
  DynFuncData* func_data = (DynFuncData*)data;

  mxArray* rhs[] = {func_data->callbacks.obj_grad, func_data->primal};

  MEX_EVAL_INTO_VEC(rhs, func_data->params, obj_grad);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
dyn_func_set_error_bound(SleqpFunc* func, double error_bound, void* data)
{
  DynFuncData* func_data = (DynFuncData*)data;

  mxArray* value = func_data->error_bound;

  *mxGetPr(value) = error_bound;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
dyn_func_set_obj_weight(SleqpFunc* func, double obj_weight, void* data)
{
  DynFuncData* func_data = (DynFuncData*)data;

  mxArray* value = func_data->obj_weight;

  *mxGetPr(value) = obj_weight;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
dyn_func_set_cons_weights(SleqpFunc* func,
                          const double* cons_weights,
                          void* data)
{
  DynFuncData* func_data = (DynFuncData*)data;

  const int num_cons = sleqp_func_num_cons(func);

  mxArray* weight_array = func_data->cons_weights;
  double* target        = mxGetPr(weight_array);

  for (int i = 0; i < num_cons; ++i)
  {
    target[i] = cons_weights[i];
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
dyn_func_cons_jac(SleqpFunc* func, SleqpSparseMatrix* cons_jac, void* data)
{
  DynFuncData* func_data = (DynFuncData*)data;

  mxArray* rhs[] = {func_data->callbacks.cons_jac, func_data->primal};

  MEX_EVAL_INTO_SPARSE_MATRIX(rhs, func_data->params, cons_jac);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
dyn_func_hess_prod(SleqpFunc* func,
                   const double* obj_dual,
                   const SleqpVec* direction,
                   const SleqpVec* cons_duals,
                   SleqpVec* product,
                   void* data)
{
  DynFuncData* func_data = (DynFuncData*)data;

  mxArray* rhs[] = {func_data->error_info};

  SLEQP_CALL(mex_hess_prod(&(func_data->hess),
                           func_data->primal,
                           obj_dual,
                           direction,
                           cons_duals,
                           rhs,
                           MEX_ARRAY_LEN(rhs),
                           product));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
dyn_func_free(void* data)
{
  DynFuncData* func_data = (DynFuncData*)data;

  SLEQP_CALL(mex_hess_free(&func_data->hess));

  mxDestroyArray(func_data->error_info);

  mxDestroyArray(func_data->cons_weights);
  mxDestroyArray(func_data->obj_weight);
  mxDestroyArray(func_data->error_bound);

  mxDestroyArray(func_data->primal);

  SLEQP_CALL(sleqp_params_release(&func_data->params));

  sleqp_free(&func_data);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
mex_dyn_func_create(SleqpFunc** star,
                    const mxArray* mex_x0,
                    const mxArray* mex_callbacks,
                    SleqpParams* params,
                    int num_vars,
                    int num_cons,
                    bool with_hess)
{
  SleqpDynFuncCallbacks callbacks = {
    .set_value        = dyn_func_set,
    .set_error_bound  = dyn_func_set_error_bound,
    .set_obj_weight   = dyn_func_set_obj_weight,
    .set_cons_weights = dyn_func_set_cons_weights,
    .eval             = dyn_func_eval,
    .obj_grad         = dyn_func_obj_grad,
    .cons_jac         = dyn_func_cons_jac,
    .hess_prod        = dyn_func_hess_prod,
    .func_free        = dyn_func_free,
  };

  DynFuncData* func_data;

  SLEQP_CALL(create_dyn_func_data(&func_data,
                                  mex_callbacks,
                                  params,
                                  num_vars,
                                  num_cons,
                                  with_hess));

  SLEQP_CALL(
    sleqp_dyn_func_create(star, &callbacks, num_vars, num_cons, func_data));

  return SLEQP_OKAY;
}
