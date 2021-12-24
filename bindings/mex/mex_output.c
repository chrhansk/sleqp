#include "mex_output.h"

#include "mex_fields.h"

static SLEQP_RETCODE
create_array_from_vec(const SleqpSparseVec* vec, mxArray** array_star)
{
  const int dim = vec->dim;

  *array_star = mxCreateDoubleMatrix(dim, 1, mxREAL);

  SLEQP_CALL(sleqp_sparse_vector_to_raw(vec, mxGetPr(*array_star)));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
set_struct_field_to_vec(mxArray* info,
                        const char* name,
                        const SleqpSparseVec* vec)
{
  mxArray* array;

  SLEQP_CALL(create_array_from_vec(vec, &array));

  mxSetField(info, 0, name, array);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
set_struct_field_to_real(mxArray* info, const char* name, double value)
{
  mxArray* array;

  array = mxCreateDoubleScalar(value);

  mxSetField(info, 0, name, array);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
create_mex_output(SleqpProblem* problem,
                  SleqpSolver* solver,
                  mxArray** sol_star,
                  mxArray** info_star)
{
  const int num_vars = sleqp_problem_num_vars(problem);

  char const* fieldnames[5] = {MEX_OUTPUT_PRIMAL,
                               MEX_OUTPUT_CONS_DUAL,
                               MEX_OUTPUT_VARS_DUAL,
                               MEX_OUTPUT_ELAPSED,
                               MEX_OUTPUT_ITER};

  *info_star = mxCreateStructMatrix(1, 1, 5, fieldnames);

  mxArray* sol  = *sol_star;
  mxArray* info = *info_star;

  SleqpIterate* iterate;

  SLEQP_CALL(sleqp_solver_solution(solver, &iterate));

  SLEQP_CALL(create_array_from_vec(sleqp_iterate_primal(iterate), sol_star));

  SLEQP_CALL(set_struct_field_to_vec(info,
                                     MEX_OUTPUT_PRIMAL,
                                     sleqp_iterate_primal(iterate)));

  SLEQP_CALL(set_struct_field_to_vec(info,
                                     MEX_OUTPUT_CONS_DUAL,
                                     sleqp_iterate_cons_dual(iterate)));

  SLEQP_CALL(set_struct_field_to_vec(info,
                                     MEX_OUTPUT_VARS_DUAL,
                                     sleqp_iterate_vars_dual(iterate)));

  SLEQP_CALL(set_struct_field_to_real(info,
                                      MEX_OUTPUT_ELAPSED,
                                      sleqp_solver_elapsed_seconds(solver)));

  SLEQP_CALL(set_struct_field_to_real(info,
                                      MEX_OUTPUT_ITER,
                                      sleqp_solver_iterations(solver)));

  return SLEQP_OKAY;
}
