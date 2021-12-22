#include <mex.h>

#include <assert.h>

#include "mex_func.h"

// TODO: Better logging
// TODO: Pass more options / params

int
array_or_cell_size(const mxArray* array)
{
  if (mxIsCell(array))
  {
    int size               = 0;
    const int num_elements = mxGetNumberOfElements(array);

    for (int i = 0; i < num_elements; ++i)
    {
      const mxArray* p               = mxGetCell(array, i);
      const int current_num_elements = mxGetNumberOfElements(p);

      size += current_num_elements;
    }

    return size;
  }
  else
  {
    return mxGetNumberOfElements(array);
  }
}

int
num_variables_from_solution(const mxArray* x0)
{
  return array_or_cell_size(x0);
}

int
num_constraints_from_options(const mxArray* options)
{
  assert(mxIsStruct(options));

  const mxArray* cons_lb = mxGetField(options, 0, "cl");

  if (cons_lb)
  {
    return array_or_cell_size(cons_lb);
  }

  return 0;
}

SLEQP_RETCODE
create_vec_from_array(SleqpSparseVec** star, const mxArray* array)
{
  const int dimension = array_or_cell_size(array);

  SLEQP_CALL(sleqp_sparse_vector_create_full(star, dimension));

  SleqpSparseVec* vec = *star;

  double* data = vec->data;

  if (mxIsCell(array))
  {
    const int num_elements = mxGetNumberOfElements(array);
    int offset             = 0;

    for (int i = 0; i < num_elements; ++i)
    {
      const mxArray* p           = mxGetCell(array, i);
      const int cur_num_elements = mxGetNumberOfElements(p);
      const double* ptr          = mxGetPr(p);

      for (int k = 0; k < cur_num_elements; ++k)
      {
        SLEQP_CALL(sleqp_sparse_vector_push(vec, k + offset, ptr[k]));
      }

      offset += cur_num_elements;
    }
  }
  else
  {
    const double* ptr = mxGetPr(array);

    for (int k = 0; k < dimension; ++k)
    {
      SLEQP_CALL(sleqp_sparse_vector_push(vec, k, ptr[k]));
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
create_vec_from_array_default(SleqpSparseVec** star,
                              const mxArray* array,
                              const int dimension,
                              double value)
{
  if (array)
  {
    SLEQP_CALL(create_vec_from_array(star, array));
  }
  else
  {
    SLEQP_CALL(sleqp_sparse_vector_create_full(star, dimension));
    SLEQP_CALL(sleqp_sparse_vector_fill(*star, value));
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
create_bounds_from_options(SleqpSparseVec** lb_star,
                           SleqpSparseVec** ub_star,
                           const int dimension,
                           const mxArray* lb_array,
                           const mxArray* ub_array)
{
  const double inf = sleqp_infinity();

  SLEQP_CALL(create_vec_from_array_default(lb_star, lb_array, dimension, -inf));

  SLEQP_CALL(create_vec_from_array_default(ub_star, ub_array, dimension, inf));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
create_var_bounds_from_options(SleqpSparseVec** var_lb_star,
                               SleqpSparseVec** var_ub_star,
                               const int num_variables,
                               const mxArray* options)
{
  return create_bounds_from_options(var_lb_star,
                                    var_ub_star,
                                    num_variables,
                                    mxGetField(options, 0, "lb"),
                                    mxGetField(options, 0, "ub"));
  return SLEQP_OKAY;
}

SLEQP_RETCODE
create_cons_bounds_from_options(SleqpSparseVec** cons_lb_star,
                                SleqpSparseVec** cons_ub_star,
                                const int num_constraints,
                                const mxArray* options)
{
  return create_bounds_from_options(cons_lb_star,
                                    cons_ub_star,
                                    num_constraints,
                                    mxGetField(options, 0, "cl"),
                                    mxGetField(options, 0, "cu"));
  return SLEQP_OKAY;
}

void
mex_log_handler(SLEQP_LOG_LEVEL level, time_t time, const char* message)
{
  mexPrintf("%s\n", message);
}

void
mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  sleqp_log_set_handler(mex_log_handler);

  if (nrhs != 3)
  {
    mexErrMsgIdAndTxt("id??",
                      "Invalid number of arguments, expected 3 found %d",
                      nrhs);
  }

  if (nlhs != 2)
  {
    mexErrMsgIdAndTxt("id??",
                      "Invalid number of return values, expected 2 found %d",
                      nlhs);
  }

  const mxArray* mex_x0      = prhs[0];
  const mxArray* mex_funcs   = prhs[1];
  const mxArray* mex_options = prhs[2];

  const int num_variables   = num_variables_from_solution(mex_x0);
  const int num_constraints = num_constraints_from_options(mex_options);

  printf("Num variables: %d, num constraints: %d\n",
         num_variables,
         num_constraints);

  SleqpSparseVec* var_lb = NULL;
  SleqpSparseVec* var_ub = NULL;

  create_var_bounds_from_options(&var_lb, &var_ub, num_variables, mex_options);

  SleqpSparseVec* cons_lb = NULL;
  SleqpSparseVec* cons_ub = NULL;

  create_cons_bounds_from_options(&cons_lb,
                                  &cons_ub,
                                  num_constraints,
                                  mex_options);

  SleqpFunc* func = NULL;

  mex_func_create(&func, mex_funcs, num_variables, num_constraints);

  SleqpProblem* problem = NULL;

  SleqpParams* params = NULL;

  sleqp_params_create(&params);

  sleqp_problem_create_simple(&problem,
                              func,
                              params,
                              var_lb,
                              var_ub,
                              cons_lb,
                              cons_ub);

  SleqpSparseVec* primal = NULL;

  create_vec_from_array(&primal, mex_x0);

  SleqpOptions* options;

  sleqp_options_create(&options);

  sleqp_options_set_int_value(options,
                              SLEQP_OPTION_INT_HESS_EVAL,
                              SLEQP_HESS_EVAL_DAMPED_BFGS);

  SleqpSolver* solver;

  sleqp_solver_create(&solver, problem, params, options, primal, NULL);

  sleqp_solver_solve(solver, SLEQP_NONE, SLEQP_NONE);

  SleqpIterate* iterate;

  sleqp_solver_solution(solver, &iterate);

  SleqpSparseVec* solution = sleqp_iterate_primal(iterate);

  plhs[0] = mxCreateDoubleMatrix(num_variables, 1, mxREAL);

  sleqp_sparse_vector_to_raw(solution, mxGetPr(plhs[0]));

  plhs[1] = mxCreateDoubleMatrix(1, 0, mxREAL);

  return;
}
