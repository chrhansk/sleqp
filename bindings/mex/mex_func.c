#include "mex_func.h"

#include <assert.h>
#include <threads.h>

#include "mex_error.h"
#include "mex_fields.h"
#include "mex_func_common.h"
#include "sleqp/pub_types.h"

static void
print_array(mxArray* array)
{
  mexCallMATLAB(0, NULL, 1, &array, MATLAB_FUNC_DISP);
}

static SLEQP_RETCODE
callback_has_field(const mxArray* mex_callbacks,
                   const char* name,
                   bool* has_field)
{
  MEX_EXPECT_STRUCT(mex_callbacks);

  mxArray* field = mxGetField(mex_callbacks, 0, name);

  (*has_field) = !!(field);

  return SLEQP_OKAY;
}

// TODO: Better handling of matlab exceptions
// TODO: Better assert

typedef struct
{
  SleqpParams* params;

  // Callbacks
  struct
  {
    mxArray* obj_val;
    mxArray* obj_grad;
    mxArray* cons_val;
    mxArray* cons_jac;
    mxArray* hess;
  } callbacks;

  bool hess_prod;

  mxArray* primal;

  mxArray* obj_dual;
  mxArray* cons_dual;
  mxArray* hess_dir;

  double* direction;
  double* product;

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

  const int nrhs = sizeof(rhs) / sizeof(rhs[0]);

  SLEQP_CALL(mex_eval_into_real(nrhs, rhs, obj_val));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
mex_func_obj_grad(SleqpFunc* func, SleqpVec* obj_grad, void* data)
{
  FuncData* func_data = (FuncData*)data;

  mxArray* rhs[] = {func_data->callbacks.obj_grad, func_data->primal};

  SLEQP_CALL(mex_eval_into_sparse_vec(2, rhs, func_data->params, obj_grad));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
mex_func_cons_val(SleqpFunc* func, SleqpVec* cons_val, void* data)
{
  FuncData* func_data = (FuncData*)data;

  mxArray* rhs[] = {func_data->callbacks.cons_val, func_data->primal};

  SLEQP_CALL(mex_eval_into_sparse_vec(2, rhs, func_data->params, cons_val));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
mex_func_cons_jac(SleqpFunc* func, SleqpSparseMatrix* cons_jac, void* data)
{
  FuncData* func_data = (FuncData*)data;

  mxArray* rhs[] = {func_data->callbacks.cons_jac, func_data->primal};

  SLEQP_CALL(mex_eval_into_sparse_matrix(2, rhs, func_data->params, cons_jac));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
hess_prod(const mxArray* hessian, const double* direction, double* product)
{
  const int dim = mxGetN(hessian);

  for (int row = 0; row < dim; ++row)
  {
    product[row] = 0.;
  }

  const mwIndex* jc = mxGetJc(hessian);
  const mwIndex* ir = mxGetIr(hessian);
  const double* pr  = mxGetPr(hessian);

  assert(jc[0] == 0);

  mwIndex index = 0;

  for (mwIndex col = 0; col < dim; ++col)
  {
    assert(jc[col] <= jc[col + 1]);

    for (; index < jc[col + 1]; ++index)
    {
      const mwIndex row  = ir[index];
      const double value = pr[index];

      assert(row >= col);

      if (row == col)
      {
        product[row] += value * direction[col];
      }
      else
      {
        product[row] += value * direction[col];
        product[col] += value * direction[row];
      }
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
hess_prod_direct(FuncData* func_data,
                 const SleqpVec* direction,
                 SleqpVec* result)
{
  SLEQP_CALL(sleqp_vec_to_raw(direction, mxGetPr(func_data->hess_dir)));

  mxArray* rhs[] = {func_data->callbacks.hess,
                    func_data->primal,
                    func_data->hess_dir,
                    func_data->obj_dual,
                    func_data->cons_dual};

  const int nrhs = sizeof(rhs) / sizeof(rhs[0]);

  SLEQP_CALL(mex_eval_into_sparse_vec(nrhs, rhs, func_data->params, result));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
hess_prod_matrix(SleqpFunc* func,
                 FuncData* func_data,
                 const SleqpVec* direction,
                 SleqpVec* result)
{
  const double zero_eps
    = sleqp_params_value(func_data->params, SLEQP_PARAM_ZERO_EPS);

  mxArray* lhs;
  mxArray* rhs[] = {func_data->callbacks.hess,
                    func_data->primal,
                    func_data->obj_dual,
                    func_data->cons_dual};

  const int nrhs = sizeof(rhs) / sizeof(rhs[0]);

  MATLAB_CALL(mexCallMATLABWithTrap(1, &lhs, nrhs, rhs, MATLAB_FUNC_FEVAL));

  MEX_EXPECT_DOUBLE(lhs);
  MEX_EXPECT_SPARSE(lhs);

  const int num_vars = sleqp_func_num_vars(func);

  MEX_EXPECT_SHAPE(lhs, num_vars, num_vars);

  SLEQP_CALL(sleqp_vec_to_raw(direction, func_data->direction));

  SLEQP_CALL(hess_prod(lhs, func_data->direction, func_data->product));

  SLEQP_CALL(
    sleqp_vec_set_from_raw(result, func_data->product, num_vars, zero_eps));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
mex_func_hess_prod(SleqpFunc* func,
                   const double* obj_dual,
                   const SleqpVec* direction,
                   const SleqpVec* cons_duals,
                   SleqpVec* result,
                   void* data)
{
  FuncData* func_data = (FuncData*)data;

  SLEQP_CALL(sleqp_vec_to_raw(cons_duals, mxGetPr(func_data->cons_dual)));

  if (obj_dual)
  {
    *mxGetPr(func_data->obj_dual) = *obj_dual;
  }
  else
  {
    *mxGetPr(func_data->obj_dual) = 0.;
  }

  if (func_data->hess_prod)
  {
    return hess_prod_direct(func_data, direction, result);
  }
  else
  {
    return hess_prod_matrix(func, func_data, direction, result);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
create_func_data(FuncData** star,
                 const mxArray* mex_callbacks,
                 SleqpParams* params,
                 int num_vars,
                 int num_cons,
                 bool with_hess)
{
  SLEQP_CALL(sleqp_malloc(star));

  FuncData* func_data = *star;

  *func_data = (FuncData){0};

  SLEQP_CALL(sleqp_params_capture(params));
  func_data->params = params;

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
    bool has_hess_prod = false;

    SLEQP_CALL(
      callback_has_field(mex_callbacks, MEX_INPUT_HESS_PROD, &has_hess_prod));

    if (has_hess_prod)
    {
      func_data->hess_prod = true;

      SLEQP_CALL(mex_callback_from_struct(mex_callbacks,
                                          MEX_INPUT_HESS_PROD,
                                          &func_data->callbacks.hess));

      func_data->hess_dir = mxCreateDoubleMatrix(num_vars, 1, mxREAL);
    }
    else
    {
      SLEQP_CALL(mex_callback_from_struct(mex_callbacks,
                                          MEX_INPUT_HESS,
                                          &func_data->callbacks.hess));

      SLEQP_CALL(sleqp_alloc_array(&func_data->direction, num_vars));
      SLEQP_CALL(sleqp_alloc_array(&func_data->product, num_vars));
    }
  }

  func_data->primal = mxCreateDoubleMatrix(1, num_vars, mxREAL);

  func_data->obj_dual  = mxCreateDoubleScalar(0.);
  func_data->cons_dual = mxCreateDoubleMatrix(num_cons, 1, mxREAL);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
mex_func_free(void* data)
{
  FuncData* func_data = (FuncData*)data;

  sleqp_free(&func_data->product);
  sleqp_free(&func_data->direction);

  mxDestroyArray(func_data->cons_dual);
  mxDestroyArray(func_data->obj_dual);
  mxDestroyArray(func_data->primal);

  if (func_data->hess_dir)
  {
    mxDestroyArray(func_data->hess_dir);
  }

  SLEQP_CALL(sleqp_params_release(&func_data->params));

  sleqp_free(&func_data);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
mex_func_create(SleqpFunc** star,
                const mxArray* mex_callbacks,
                SleqpParams* params,
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
                              params,
                              num_vars,
                              num_cons,
                              with_hess));

  SLEQP_CALL(
    sleqp_func_create(star, &callbacks, num_vars, num_cons, (void*)func_data));

  return SLEQP_OKAY;
}
