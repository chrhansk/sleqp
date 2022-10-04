#include "mex_func_common.h"

#include "mex_error.h"
#include "mex_fields.h"

SLEQP_RETCODE
mex_callback_from_struct(const mxArray* mex_callbacks,
                         const char* name,
                         mxArray** star)
{
  MEX_EXPECT_STRUCT(mex_callbacks);

  *star = mxGetField(mex_callbacks, 0, name);

  MEX_EXPECT_FUNCTION_HANDLE(*star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
mex_callback_has_field(const mxArray* mex_callbacks,
                       const char* name,
                       bool* has_field)
{
  MEX_EXPECT_STRUCT(mex_callbacks);

  mxArray* field = mxGetField(mex_callbacks, 0, name);

  (*has_field) = !!(field);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
mex_eval_into_real(int nrhs, mxArray** rhs, double* value)
{
  mxArray* lhs[] = {NULL};

  MEX_CALL(mexCallMATLABWithTrap(MEX_ARRAY_LEN(lhs),
                                 lhs,
                                 nrhs,
                                 rhs,
                                 MATLAB_FUNC_FEVAL));

  SLEQP_CALL(mex_array_to_real(*lhs, value));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
mex_array_to_real(mxArray* array, double* value)
{
  MEX_EXPECT_DOUBLE(array);
  MEX_EXPECT_SCALAR(array);

  *value = *mxGetPr(array);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
mex_eval_into_bool(int nrhs, mxArray** rhs, bool* value)
{
  mxArray* lhs[] = {NULL};

  MEX_CALL(mexCallMATLABWithTrap(MEX_ARRAY_LEN(lhs),
                                 lhs,
                                 nrhs,
                                 rhs,
                                 MATLAB_FUNC_FEVAL));

  MEX_EXPECT_LOGICAL_SCALAR(*lhs);

  *value = mxIsLogicalScalarTrue(*lhs);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
mex_eval_into_vec(int nrhs, mxArray** rhs, SleqpParams* params, SleqpVec* vec)
{
  mxArray* lhs[] = {NULL};

  MEX_CALL(mexCallMATLABWithTrap(MEX_ARRAY_LEN(lhs),
                                 lhs,
                                 nrhs,
                                 rhs,
                                 MATLAB_FUNC_FEVAL));

  SLEQP_CALL(mex_array_to_vec(*lhs, params, vec));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
mex_array_to_vec(const mxArray* array, SleqpParams* params, SleqpVec* vec)
{
  const double zero_eps = sleqp_params_value(params, SLEQP_PARAM_ZERO_EPS);

  MEX_EXPECT_DOUBLE(array);
  MEX_EXPECT_NUM_ELEMENTS(array, vec->dim);

  SLEQP_CALL(sleqp_vec_set_from_raw(vec, mxGetPr(array), vec->dim, zero_eps));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
array_to_sparse_matrix(const mxArray* array, SleqpSparseMatrix* matrix)
{
  MEX_EXPECT_DOUBLE(array);
  MEX_EXPECT_SPARSE(array);

  const int num_cols = sleqp_sparse_matrix_num_cols(matrix);
  const int num_rows = sleqp_sparse_matrix_num_rows(matrix);

  MEX_EXPECT_SHAPE(array, num_rows, num_cols);

  const mwIndex* jc = mxGetJc(array);
  const mwIndex* ir = mxGetIr(array);
  const double* pr  = mxGetPr(array);

  MEX_EXPECT_NOT_NULL(jc);
  MEX_EXPECT_NOT_NULL(ir);
  MEX_EXPECT_NOT_NULL(pr);

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
  mxArray* lhs[] = {NULL};

  MEX_CALL(mexCallMATLABWithTrap(MEX_ARRAY_LEN(lhs),
                                 lhs,
                                 nrhs,
                                 rhs,
                                 MATLAB_FUNC_FEVAL));

  SLEQP_CALL(array_to_sparse_matrix(*lhs, matrix));

  return SLEQP_OKAY;
}
