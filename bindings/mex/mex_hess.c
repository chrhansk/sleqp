#include "mex_hess.h"

#include "mex_error.h"
#include "mex_fields.h"
#include "mex_func_common.h"

SLEQP_RETCODE
mex_hess_init(MexHess* hess,
              SleqpSettings* settings,
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

  SLEQP_CALL(sleqp_settings_capture(settings));
  hess->settings = settings;

  hess->cons_dual = mxCreateDoubleMatrix(num_cons, 1, mxREAL);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
reserve_args(MexHess* hess, int nargs)
{
  if (hess->nargs > nargs)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_realloc(&hess->args, nargs));
  hess->nargs = nargs;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
prepare_args(MexHess* hess,
             mxArray** default_args,
             int ndefault_args,
             mxArray** extra_args,
             int nextra_args)
{
  const int ntotal_args = ndefault_args + nextra_args;

  reserve_args(hess, ntotal_args);

  {
    for (int i = 0; i < ndefault_args; ++i)
    {
      hess->args[i] = default_args[i];
    }

    for (int i = 0; i < nextra_args; ++i)
    {
      hess->args[ndefault_args + i] = extra_args[i];
    }
  }

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

  MEX_EXPECT_NOT_NULL(jc);
  MEX_EXPECT_NOT_NULL(ir);
  MEX_EXPECT_NOT_NULL(pr);

  assert(jc[0] == 0);

  mwIndex index = 0;

  for (mwIndex col = 0; col < dim; ++col)
  {
    assert(jc[col] <= jc[col + 1]);

    for (; index < jc[col + 1]; ++index)
    {
      const mwIndex row  = ir[index];
      const double value = pr[index];

      if (row < col)
      {
        sleqp_raise(SLEQP_FUNC_EVAL_ERROR,
                    "Hessian entry at (%ld, %ld) above diagonal",
                    row,
                    col);
      }

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
                 mxArray** rhs,
                 int nrhs,
                 SleqpVec* result)
{
  SLEQP_CALL(sleqp_vec_to_raw(direction, mxGetPr(hess->hess_dir)));

  mxArray* default_args[] = {hess->callbacks.hess,
                             primal,
                             hess->hess_dir,
                             hess->cons_dual};

  if (nrhs == 0)
  {
    MEX_EVAL_INTO_VEC(default_args, hess->settings, result);
  }
  else
  {
    const int ndefault_args = MEX_ARRAY_LEN(default_args);
    const int ntotal_args   = ndefault_args + nrhs;

    SLEQP_CALL(prepare_args(hess, default_args, ndefault_args, rhs, nrhs));

    SLEQP_CALL(
      mex_eval_into_vec(ntotal_args, hess->args, hess->settings, result));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
hess_prod_matrix(MexHess* hess,
                 mxArray* primal,
                 const SleqpVec* direction,
                 mxArray** rhs,
                 int nrhs,
                 SleqpVec* result)
{
  const double zero_eps
    = sleqp_settings_real_value(hess->settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  mxArray* lhs[] = {NULL};

  mxArray* default_args[]
    = {hess->callbacks.hess, primal, hess->cons_dual};

  if (nrhs == 0)
  {
    MEX_EVAL(lhs, default_args);
  }
  else
  {
    const int ndefault_args = MEX_ARRAY_LEN(default_args);
    const int ntotal_args   = ndefault_args + nrhs;

    SLEQP_CALL(prepare_args(hess, default_args, ndefault_args, rhs, nrhs));

    MEX_CALL(mexCallMATLABWithTrap(MEX_ARRAY_LEN(lhs),
                                   lhs,
                                   ntotal_args,
                                   hess->args,
                                   MATLAB_FUNC_FEVAL));
  }

  MEX_EXPECT_DOUBLE(lhs[0]);
  MEX_EXPECT_SPARSE(lhs[0]);

  const int num_vars = mxGetNumberOfElements(primal);

  MEX_EXPECT_SHAPE(lhs[0], num_vars, num_vars);

  SLEQP_CALL(sleqp_vec_to_raw(direction, hess->direction));

  SLEQP_CALL(prod_from_hess_matrix(lhs[0], hess->direction, hess->product));

  SLEQP_CALL(sleqp_vec_set_from_raw(result, hess->product, num_vars, zero_eps));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
mex_hess_prod(MexHess* hess,
              mxArray* primal,
              const SleqpVec* direction,
              const SleqpVec* cons_duals,
              mxArray** rhs,
              int nrhs,
              SleqpVec* result)
{
  SLEQP_CALL(sleqp_vec_to_raw(cons_duals, mxGetPr(hess->cons_dual)));

  if (!!(hess->hess_dir))
  {
    return hess_prod_direct(hess, primal, direction, rhs, nrhs, result);
  }
  else
  {
    return hess_prod_matrix(hess, primal, direction, rhs, nrhs, result);
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
mex_hess_free(MexHess* hess)
{
  sleqp_free(&hess->args);

  mxDestroyArray(hess->cons_dual);

  sleqp_free(&hess->product);
  sleqp_free(&hess->direction);

  if (hess->hess_dir)
  {
    mxDestroyArray(hess->hess_dir);
  }

  SLEQP_CALL(sleqp_settings_release(&hess->settings));

  return SLEQP_OKAY;
}
