#include "mex_func.h"

#include <assert.h>

#define MEX_CALL(x)                                                            \
  do                                                                           \
  {                                                                            \
    mxArray* exception = (x);                                                  \
                                                                               \
    if (exception)                                                             \
    {                                                                          \
      return SLEQP_INTERNAL_ERROR;                                             \
    }                                                                          \
                                                                               \
  } while (0)

// TODO: Better handling of matlab exceptions
// TODO: Perform dimension checks
// TODO: Better assert
// TODO: More type checks
// TODO: Add zero eps

typedef struct
{
  // Callbacks
  mxArray* objective;
  mxArray* gradient;
  mxArray* constraints;
  mxArray* jacobian;
  mxArray* hessian;

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
  mxArray* rhs[] = {func_data->objective, func_data->primal};

  MEX_CALL(mexCallMATLABWithTrap(1, &lhs, 2, rhs, "feval"));

  assert(mxIsDouble(lhs));

  *obj_val = *mxGetPr(lhs);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
mex_func_obj_grad(SleqpFunc* func, SleqpSparseVec* obj_grad, void* data)
{
  FuncData* func_data = (FuncData*)data;

  mxArray* lhs;
  mxArray* rhs[] = {func_data->gradient, func_data->primal};

  MEX_CALL(mexCallMATLABWithTrap(1, &lhs, 2, rhs, "feval"));

  const int num_vars = sleqp_func_num_vars(func);

  assert(mxIsDouble(lhs));
  assert(mxGetNumberOfElements(lhs) == num_vars);

  SLEQP_CALL(
    sleqp_sparse_vector_from_raw(obj_grad, mxGetPr(lhs), num_vars, 0.));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
mex_func_cons_val(SleqpFunc* func,
                  const SleqpSparseVec* cons_indices,
                  SleqpSparseVec* cons_val,
                  void* data)
{
  FuncData* func_data = (FuncData*)data;

  mxArray* lhs;
  mxArray* rhs[] = {func_data->constraints, func_data->primal};

  MEX_CALL(mexCallMATLABWithTrap(1, &lhs, 2, rhs, "feval"));

  const int num_cons = sleqp_func_num_cons(func);

  assert(mxGetNumberOfElements(lhs) == num_cons);

  SLEQP_CALL(
    sleqp_sparse_vector_from_raw(cons_val, mxGetPr(lhs), num_cons, 0.));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
array_to_sparse_matrix(const mxArray* array, SleqpSparseMatrix* matrix)
{
  assert(mxIsSparse(array));

  const int num_cols = sleqp_sparse_matrix_num_cols(matrix);
  const int num_rows = sleqp_sparse_matrix_num_rows(matrix);

  assert(mxGetM(array) == num_rows);
  assert(mxGetN(array) == num_cols);

  const mwIndex* jc = mxGetJc(array);
  const mwIndex* ir = mxGetIr(array);
  const double* pr  = mxGetPr(array);

  const int nnz = jc[num_cols];

  SLEQP_CALL(sleqp_sparse_matrix_reserve(matrix, nnz));

  assert(jc[0] == 0);

  mwIndex index = 0;

  for (mwIndex col = 0; col < num_cols; ++col)
  {
    SLEQP_CALL(sleqp_sparse_matrix_push_column(matrix, col));

    assert(jc[col] <= jc[col + 1]);

    for (; index < jc[col + 1]; ++index)
    {
      const mwIndex row  = ir[index];
      const double value = pr[index];

      SLEQP_CALL(sleqp_sparse_matrix_push(matrix, row, col, value));
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
mex_func_cons_jac(SleqpFunc* func,
                  const SleqpSparseVec* cons_indices,
                  SleqpSparseMatrix* cons_jac,
                  void* data)
{
  FuncData* func_data = (FuncData*)data;

  mxArray* lhs;
  mxArray* rhs[] = {func_data->jacobian, func_data->primal};

  MEX_CALL(mexCallMATLABWithTrap(1, &lhs, 2, rhs, "feval"));

  SLEQP_CALL(array_to_sparse_matrix(lhs, cons_jac));

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
  mxArray* rhs[] = {func_data->hessian,
                    func_data->primal,
                    func_data->obj_dual,
                    func_data->cons_dual};

  MEX_CALL(mexCallMATLABWithTrap(1, &lhs, 4, rhs, "feval"));

  assert(mxIsDouble(lhs));
  assert(!mxIsComplex(lhs));
  assert(mxIsSparse(lhs));

  const int num_vars = sleqp_func_num_vars(func);

  assert(mxGetM(lhs) == num_vars);
  assert(mxGetN(lhs) == num_vars);

  SLEQP_CALL(sleqp_sparse_vector_to_raw(direction, func_data->direction));

  SLEQP_CALL(hess_prod(lhs, func_data->direction, func_data->product));

  SLEQP_CALL(
    sleqp_sparse_vector_from_raw(result, func_data->product, num_vars, 0.));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
create_func_data(FuncData** star,
                 const mxArray* callbacks,
                 int num_vars,
                 int num_cons)
{
  SLEQP_CALL(sleqp_malloc(star));

  FuncData* func_data = *star;

  *func_data = (FuncData){0};

  assert(mxIsStruct(callbacks));

  func_data->objective   = mxGetField(callbacks, 0, "objective");
  func_data->gradient    = mxGetField(callbacks, 0, "gradient");
  func_data->constraints = mxGetField(callbacks, 0, "constraints");
  func_data->jacobian    = mxGetField(callbacks, 0, "jacobian");
  func_data->hessian     = mxGetField(callbacks, 0, "hessian");

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

  sleqp_free(&func_data);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
mex_func_create(SleqpFunc** star,
                const mxArray* mex_callbacks,
                int num_variables,
                int num_constraints)
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
                              num_variables,
                              num_constraints));

  SLEQP_CALL(sleqp_func_create(star,
                               &callbacks,
                               num_variables,
                               num_constraints,
                               (void*)func_data));

  return SLEQP_OKAY;
}
