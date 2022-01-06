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

static SLEQP_RETCODE
create_working_set_output(mxArray* info,
                          SleqpProblem* problem,
                          const SleqpWorkingSet* working_set)
{
  const int num_vars = sleqp_problem_num_vars(problem);
  const int num_cons = sleqp_problem_num_cons(problem);

  {
    mxArray* array = mxCreateDoubleMatrix(num_vars, 1, mxREAL);

    double* values = mxGetPr(array);

    for (int j = 0; j < num_vars; ++j)
    {
      values[j] = sleqp_working_set_var_state(working_set, j);
    }

    mxSetField(info, 0, MEX_OUTPUT_WORKING_VARS, array);
  }

  {
    mxArray* array = mxCreateDoubleMatrix(num_cons, 1, mxREAL);

    double* values = mxGetPr(array);

    for (int i = 0; i < num_cons; ++i)
    {
      values[i] = sleqp_working_set_cons_state(working_set, i);
    }

    mxSetField(info, 0, MEX_OUTPUT_WORKING_CONS, array);
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
mex_create_solver_output(SleqpProblem* problem,
                         SleqpSolver* solver,
                         mxArray** sol_star,
                         mxArray** info_star)
{
  const char* fieldnames[] = {MEX_OUTPUT_PRIMAL,
                              MEX_OUTPUT_CONS_DUAL,
                              MEX_OUTPUT_VARS_DUAL,
                              MEX_OUTPUT_ELAPSED,
                              MEX_OUTPUT_ITER,
                              MEX_OUTPUT_STATUS,
                              MEX_OUTPUT_WORKING_VARS,
                              MEX_OUTPUT_WORKING_CONS};

  const int num_fields = sizeof(fieldnames) / sizeof(const char*);

  *info_star = mxCreateStructMatrix(1, 1, num_fields, fieldnames);

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

  SLEQP_CALL(set_struct_field_to_real(info,
                                      MEX_OUTPUT_STATUS,
                                      sleqp_solver_status(solver)));

  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);

  SLEQP_CALL(create_working_set_output(info, problem, working_set));

  return SLEQP_OKAY;
}
