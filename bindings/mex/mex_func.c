#include "mex_func.h"

#include <assert.h>
#include <threads.h>

#include "mex_fields.h"
#include "mex_func_common.h"

static void
print_array(mxArray* array)
{
  mexCallMATLAB(0, NULL, 1, &array, MATLAB_FUNC_DISP);
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
    mxArray* hessian;
  } callbacks;

  mxArray* primal;

  mxArray* obj_dual;
  mxArray* cons_dual;

  double* direction;
  double* product;

} FuncData;

static SLEQP_RETCODE
mex_func_set(SleqpFunc* func,
             SleqpSparseVec* x,
             SLEQP_VALUE_REASON reason,
             bool* reject,
             int* obj_grad_nnz,
             int* cons_val_nnz,
             int* cons_jac_nnz,
             void* data)
{
  FuncData* func_data = (FuncData*)data;

  *obj_grad_nnz = 0;
  *cons_val_nnz = 0;
  *cons_jac_nnz = 0;

  SLEQP_CALL(sleqp_sparse_vector_to_raw(x, mxGetPr(func_data->primal)));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
mex_func_obj_val(SleqpFunc* func, double* obj_val, void* data)
{
  FuncData* func_data = (FuncData*)data;

  mxArray* lhs;
  mxArray* rhs[] = {func_data->callbacks.obj_val, func_data->primal};

  MATLAB_CALL(mexCallMATLABWithTrap(1, &lhs, 2, rhs, MATLAB_FUNC_FEVAL));

  assert(mxIsScalar(lhs));
  assert(!mxIsComplex(lhs));
  assert(mxIsDouble(lhs));

  *obj_val = *mxGetPr(lhs);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
mex_func_obj_grad(SleqpFunc* func, SleqpSparseVec* obj_grad, void* data)
{
  FuncData* func_data = (FuncData*)data;

  mxArray* rhs[] = {func_data->callbacks.obj_grad, func_data->primal};

  SLEQP_CALL(mex_eval_into_sparse_vec(2, rhs, func_data->params, obj_grad));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
mex_func_cons_val(SleqpFunc* func,
                  const SleqpSparseVec* cons_indices,
                  SleqpSparseVec* cons_val,
                  void* data)
{
  FuncData* func_data = (FuncData*)data;

  mxArray* rhs[] = {func_data->callbacks.cons_val, func_data->primal};

  SLEQP_CALL(mex_eval_into_sparse_vec(2, rhs, func_data->params, cons_val));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
mex_func_cons_jac(SleqpFunc* func,
                  const SleqpSparseVec* cons_indices,
                  SleqpSparseMatrix* cons_jac,
                  void* data)
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
mex_func_hess_prod(SleqpFunc* func,
                   const double* obj_dual,
                   const SleqpSparseVec* direction,
                   const SleqpSparseVec* cons_duals,
                   SleqpSparseVec* result,
                   void* data)
{
  FuncData* func_data = (FuncData*)data;

  const double zero_eps
    = sleqp_params_value(func_data->params, SLEQP_PARAM_ZERO_EPS);

  SLEQP_CALL(
    sleqp_sparse_vector_to_raw(cons_duals, mxGetPr(func_data->cons_dual)));

  if (obj_dual)
  {
    *mxGetPr(func_data->obj_dual) = *obj_dual;
  }
  else
  {
    *mxGetPr(func_data->obj_dual) = 0.;
  }

  mxArray* lhs;
  mxArray* rhs[] = {func_data->callbacks.hessian,
                    func_data->primal,
                    func_data->obj_dual,
                    func_data->cons_dual};

  MATLAB_CALL(mexCallMATLABWithTrap(1, &lhs, 4, rhs, MATLAB_FUNC_FEVAL));

  assert(mxIsDouble(lhs));
  assert(!mxIsComplex(lhs));
  assert(mxIsSparse(lhs));

  const int num_vars = sleqp_func_num_vars(func);

  assert(mxGetM(lhs) == num_vars);
  assert(mxGetN(lhs) == num_vars);

  SLEQP_CALL(sleqp_sparse_vector_to_raw(direction, func_data->direction));

  SLEQP_CALL(hess_prod(lhs, func_data->direction, func_data->product));

  SLEQP_CALL(sleqp_sparse_vector_from_raw(result,
                                          func_data->product,
                                          num_vars,
                                          zero_eps));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
create_func_data(FuncData** star,
                 const mxArray* mex_callbacks,
                 bool with_hessian,
                 SleqpParams* params,
                 int num_vars,
                 int num_cons)
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

  if (with_hessian)
  {
    SLEQP_CALL(mex_callback_from_struct(mex_callbacks,
                                        MEX_INPUT_HESS,
                                        &func_data->callbacks.hessian));
  }

  func_data->primal = mxCreateDoubleMatrix(1, num_vars, mxREAL);

  func_data->obj_dual  = mxCreateDoubleScalar(0.);
  func_data->cons_dual = mxCreateDoubleMatrix(num_cons, 1, mxREAL);

  SLEQP_CALL(sleqp_alloc_array(&func_data->direction, num_vars));
  SLEQP_CALL(sleqp_alloc_array(&func_data->product, num_vars));

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

  SLEQP_CALL(sleqp_params_release(&func_data->params));

  sleqp_free(&func_data);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
mex_func_create(SleqpFunc** star,
                const mxArray* mex_callbacks,
                bool with_hessian,
                SleqpParams* params,
                int num_variables,
                int num_constraints)
{
  SleqpFuncCallbacks callbacks
    = {.set_value = mex_func_set,
       .obj_val   = mex_func_obj_val,
       .obj_grad  = mex_func_obj_grad,
       .cons_val  = mex_func_cons_val,
       .cons_jac  = mex_func_cons_jac,
       .hess_prod = with_hessian ? mex_func_hess_prod : NULL,
       .func_free = mex_func_free};

  FuncData* func_data;

  SLEQP_CALL(create_func_data(&func_data,
                              mex_callbacks,
                              with_hessian,
                              params,
                              num_variables,
                              num_constraints));

  SLEQP_CALL(sleqp_func_create(star,
                               &callbacks,
                               num_variables,
                               num_constraints,
                               (void*)func_data));

  return SLEQP_OKAY;
}
