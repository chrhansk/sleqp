#include "mex_solve_common.h"

#include <assert.h>
#include <mex.h>

#include "mex_fields.h"
#include "mex_output.h"
#include "mex_problem.h"

typedef struct
{
  const char* name;
  SLEQP_PARAM param;
} ParamName;

static const ParamName param_names[] = {
  {MEX_PARAM_ZERO_EPS, SLEQP_PARAM_ZERO_EPS},
  {MEX_PARAM_EPS, SLEQP_PARAM_EPS},
  {MEX_PARAM_OBJ_LOWER, SLEQP_PARAM_OBJ_LOWER},
  {MEX_PARAM_DERIV_PERTURBATION, SLEQP_PARAM_DERIV_PERTURBATION},
  {MEX_PARAM_DERIV_TOL, SLEQP_PARAM_DERIV_TOL},
  {MEX_PARAM_CAUCHY_TAU, SLEQP_PARAM_CAUCHY_TAU},
  {MEX_PARAM_CAUCHY_ETA, SLEQP_PARAM_CAUCHY_ETA},
  {MEX_PARAM_LINESEARCH_TAU, SLEQP_PARAM_LINESEARCH_TAU},
  {MEX_PARAM_LINESEARCH_ETA, SLEQP_PARAM_LINESEARCH_ETA},
  {MEX_PARAM_LINESEARCH_CUTOFF, SLEQP_PARAM_LINESEARCH_CUTOFF},
  {MEX_PARAM_FEASIBILITY_TOL, SLEQP_PARAM_FEASIBILITY_TOL},
  {MEX_PARAM_SLACKNESS_TOL, SLEQP_PARAM_SLACKNESS_TOL},
  {MEX_PARAM_STATIONARITY_TOL, SLEQP_PARAM_STATIONARITY_TOL},
  {MEX_PARAM_ACCEPTED_REDUCTION, SLEQP_PARAM_ACCEPTED_REDUCTION},
  {MEX_PARAM_DEADPOINT_BOUND, SLEQP_PARAM_DEADPOINT_BOUND},
};

static SLEQP_RETCODE
read_params(SleqpParams* params, const mxArray* mex_options)
{
  if (!mxIsStruct(mex_options))
  {
    return SLEQP_ILLEGAL_ARGUMENT;
  }

  const int num_params = sizeof(param_names) / sizeof(param_names[0]);

  for (int i = 0; i < num_params; ++i)
  {
    const ParamName* param_name = param_names + i;
    const mxArray* value        = mxGetField(mex_options, 0, param_name->name);

    if (!value)
    {
      continue;
    }

    if (!(mxIsScalar(value) && mxIsDouble(value)))
    {
      return SLEQP_ILLEGAL_ARGUMENT;
    }

    const double* param_ptr = mxGetPr(value);

    assert(param_ptr);

    const double param_value = *param_ptr;

    SLEQP_CALL(sleqp_params_set_value(params, param_name->param, param_value));
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
mex_solve(mxArray** sol_star,
          mxArray** info_star,
          bool lsq,
          const mxArray* mex_x0,
          const mxArray* mex_funcs,
          const mxArray* mex_options)
{
  SleqpOptions* options;
  SleqpParams* params;
  SleqpProblem* problem;
  SleqpSolver* solver;
  SleqpSparseVec* initial;

  SLEQP_CALL(sleqp_options_create(&options));
  SLEQP_CALL(sleqp_params_create(&params));

  SLEQP_CALL(read_params(params, mex_options));

  SLEQP_CALL(
    mex_problem_create(&problem, params, lsq, mex_x0, mex_funcs, mex_options));

  SLEQP_CALL(mex_create_vec_from_array(&initial, mex_x0));

  SLEQP_CALL(
    sleqp_solver_create(&solver, problem, params, options, initial, NULL));

  SLEQP_CALL(sleqp_solver_solve(solver, SLEQP_NONE, SLEQP_NONE));

  SLEQP_CALL(mex_create_solver_output(problem, solver, sol_star, info_star));

  SLEQP_CALL(sleqp_sparse_vector_free(&initial));
  SLEQP_CALL(sleqp_solver_release(&solver));
  SLEQP_CALL(sleqp_problem_release(&problem));

  SLEQP_CALL(sleqp_params_release(&params));
  SLEQP_CALL(sleqp_options_release(&options));

  return SLEQP_OKAY;
}
