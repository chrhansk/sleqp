#include "mex_func_common.h"

#include "mex_fields.h"

SLEQP_RETCODE
mex_callback_from_struct(const mxArray* mex_callbacks,
                         const char* name,
                         mxArray** star)
{
  if (!mxIsStruct(mex_callbacks))
  {
    return SLEQP_ERROR;
  }

  *star = mxGetField(mex_callbacks, 0, name);

  if (!(*star && mxIsFunctionHandle(*star)))
  {
    return SLEQP_ERROR;
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
mex_callback_has_field(const mxArray* mex_callbacks,
                       const char* name,
                       bool* has_field)
{
  if (!mxIsStruct(mex_callbacks))
  {
    return SLEQP_ERROR;
  }

  mxArray* field = mxGetField(mex_callbacks, 0, name);

  (*has_field) = !!(field);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
mex_eval_into_real(int nrhs, mxArray** rhs, double* value)
{
  mxArray* lhs;

  MATLAB_CALL(mexCallMATLABWithTrap(1, &lhs, nrhs, rhs, MATLAB_FUNC_FEVAL));

  if (!mxIsDouble(lhs) || mxIsComplex(lhs))
  {
    return SLEQP_ERROR;
  }

  if (!mxIsScalar(lhs))
  {
    return SLEQP_ERROR;
  }

  *value = *mxGetPr(lhs);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
mex_eval_into_bool(int nrhs, mxArray** rhs, bool* value)
{
  mxArray* lhs;

  MATLAB_CALL(mexCallMATLABWithTrap(1, &lhs, nrhs, rhs, MATLAB_FUNC_FEVAL));

  if (!mxIsLogicalScalar(lhs))
  {
    return SLEQP_ERROR;
  }

  *value = mxIsLogicalScalarTrue(lhs);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
mex_eval_into_sparse_vec(int nrhs,
                         mxArray** rhs,
                         SleqpParams* params,
                         SleqpSparseVec* vec)
{
  const double zero_eps = sleqp_params_value(params, SLEQP_PARAM_ZERO_EPS);

  mxArray* lhs;

  MATLAB_CALL(mexCallMATLABWithTrap(1, &lhs, nrhs, rhs, MATLAB_FUNC_FEVAL));

  if (!mxIsDouble(lhs) || mxIsComplex(lhs))
  {
    return SLEQP_ERROR;
  }

  if (mxGetNumberOfElements(lhs) != vec->dim)
  {
    return SLEQP_ERROR;
  }

  SLEQP_CALL(
    sleqp_sparse_vector_from_raw(vec, mxGetPr(lhs), vec->dim, zero_eps));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
array_to_sparse_matrix(const mxArray* array, SleqpSparseMatrix* matrix)
{
  assert(mxIsSparse(array));
  assert(mxIsDouble(array));
  assert(!mxIsComplex(array));

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

SLEQP_RETCODE
mex_eval_into_sparse_matrix(int nrhs,
                            mxArray** rhs,
                            SleqpParams* params,
                            SleqpSparseMatrix* matrix)
{
  mxArray* lhs;

  MATLAB_CALL(mexCallMATLABWithTrap(1, &lhs, nrhs, rhs, MATLAB_FUNC_FEVAL));

  if (!mxIsDouble(lhs) || mxIsComplex(lhs) || !mxIsSparse(lhs))
  {
    return SLEQP_ERROR;
  }

  SLEQP_CALL(array_to_sparse_matrix(lhs, matrix));

  return SLEQP_OKAY;
}
