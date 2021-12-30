#include "mex_solve_common.h"

#include "mex_output.h"
#include "mex_problem.h"

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
