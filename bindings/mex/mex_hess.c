#include "mex_hess.h"

#include "mex_error.h"
#include "mex_fields.h"
#include "mex_func_common.h"

SLEQP_RETCODE
mex_hess_init(MexHess* hess,
              SleqpParams* params,
              const mxArray* mex_callbacks,
              const int num_vars,
              const int num_cons)
{
  bool has_hess_prod = false;

  SLEQP_CALL(
    mex_callback_has_field(mex_callbacks, MEX_INPUT_HESS_PROD, &has_hess_prod));

  if (has_hess_prod)
  {
    SLEQP_CALL(mex_callback_from_struct(mex_callbacks,
                                        MEX_INPUT_HESS_PROD,
                                        &hess->callbacks.hess));

    hess->hess_dir = mxCreateDoubleMatrix(num_vars, 1, mxREAL);
  }
  else
  {
    SLEQP_CALL(mex_callback_from_struct(mex_callbacks,
                                        MEX_INPUT_HESS,
                                        &hess->callbacks.hess));

    SLEQP_CALL(sleqp_alloc_array(&hess->direction, num_vars));
    SLEQP_CALL(sleqp_alloc_array(&hess->product, num_vars));
  }

  SLEQP_CALL(sleqp_params_capture(params));
  hess->params = params;

  hess->obj_dual  = mxCreateDoubleScalar(0.);
  hess->cons_dual = mxCreateDoubleMatrix(num_cons, 1, mxREAL);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
prod_from_hess_matrix(const mxArray* hessian,
                      const double* direction,
                      double* product)
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
hess_prod_direct(MexHess* hess,
                 mxArray* primal,
                 const SleqpVec* direction,
                 SleqpVec* result)
{
  SLEQP_CALL(sleqp_vec_to_raw(direction, mxGetPr(hess->hess_dir)));

  mxArray* rhs[] = {hess->callbacks.hess,
                    primal,
                    hess->hess_dir,
                    hess->obj_dual,
                    hess->cons_dual};

  const int nrhs = sizeof(rhs) / sizeof(rhs[0]);

  SLEQP_CALL(mex_eval_into_sparse_vec(nrhs, rhs, hess->params, result));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
hess_prod_matrix(MexHess* hess,
                 mxArray* primal,
                 const SleqpVec* direction,
                 SleqpVec* result)
{
  const double zero_eps
    = sleqp_params_value(hess->params, SLEQP_PARAM_ZERO_EPS);

  mxArray* lhs;
  mxArray* rhs[]
    = {hess->callbacks.hess, primal, hess->obj_dual, hess->cons_dual};

  const int nrhs = sizeof(rhs) / sizeof(rhs[0]);

  MATLAB_CALL(mexCallMATLABWithTrap(1, &lhs, nrhs, rhs, MATLAB_FUNC_FEVAL));

  MEX_EXPECT_DOUBLE(lhs);
  MEX_EXPECT_SPARSE(lhs);

  const int num_vars = mxGetNumberOfElements(primal);

  MEX_EXPECT_SHAPE(lhs, num_vars, num_vars);

  SLEQP_CALL(sleqp_vec_to_raw(direction, hess->direction));

  SLEQP_CALL(prod_from_hess_matrix(lhs, hess->direction, hess->product));

  SLEQP_CALL(sleqp_vec_set_from_raw(result, hess->product, num_vars, zero_eps));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
mex_hess_prod(MexHess* hess,
              mxArray* primal,
              const double* obj_dual,
              const SleqpVec* direction,
              const SleqpVec* cons_duals,
              SleqpVec* result)
{
  SLEQP_CALL(sleqp_vec_to_raw(cons_duals, mxGetPr(hess->cons_dual)));

  if (obj_dual)
  {
    *mxGetPr(hess->obj_dual) = *obj_dual;
  }
  else
  {
    *mxGetPr(hess->obj_dual) = 0.;
  }

  if (!!(hess->hess_dir))
  {
    return hess_prod_direct(hess, primal, direction, result);
  }
  else
  {
    return hess_prod_matrix(hess, primal, direction, result);
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
mex_hess_free(MexHess* hess)
{
  mxDestroyArray(hess->cons_dual);
  mxDestroyArray(hess->obj_dual);

  sleqp_free(&hess->product);
  sleqp_free(&hess->direction);

  if (hess->hess_dir)
  {
    mxDestroyArray(hess->hess_dir);
  }

  SLEQP_CALL(sleqp_params_release(&hess->params));

  return SLEQP_OKAY;
}
