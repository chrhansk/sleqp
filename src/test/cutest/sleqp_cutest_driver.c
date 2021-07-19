#include "sleqp_cutest_driver.h"

#include <assert.h>

#include "iterate.h"
#include "log.h"
#include "mem.h"
#include "options.h"

#include "sleqp_cutest_defs.h"
#include "sleqp_cutest_data.h"
#include "sleqp_cutest_constrained.h"
#include "sleqp_cutest_unconstrained.h"

static
SLEQP_RETCODE report_result(SleqpSolver* solver,
                            SleqpProblem* problem,
                            const char* probname)
{
  const char* descriptions[] = {
    [SLEQP_FEASIBLE] = "feasible",
    [SLEQP_OPTIMAL] = "optimal",
    [SLEQP_INFEASIBLE] = "infeasible",
    [SLEQP_INVALID] = "invalid"
  };

  const int num_variables = sleqp_problem_num_variables(problem);
  const int num_constraints = sleqp_problem_num_constraints(problem);

  SLEQP_STATUS status = sleqp_solver_get_status(solver);
  SleqpIterate* iterate;

  SLEQP_CALL(sleqp_solver_get_solution(solver, &iterate));

  const int iterations = sleqp_solver_get_iterations(solver);

  int last_step_bdry;

  SLEQP_CALL(sleqp_solver_get_int_state(solver,
                                        SLEQP_SOLVER_STATE_INT_LAST_STEP_ON_BDRY,
                                        &last_step_bdry));

  double last_trust_radius;

  SLEQP_CALL(sleqp_solver_get_real_state(solver,
                                         SLEQP_SOLVER_STATE_REAL_TRUST_RADIUS,
                                         &last_trust_radius));

  double min_rayleigh;

  SLEQP_CALL(sleqp_solver_get_real_state(solver,
                                         SLEQP_SOLVER_STATE_REAL_MIN_RAYLEIGH,
                                         &min_rayleigh));

  double max_rayleigh;

  SLEQP_CALL(sleqp_solver_get_real_state(solver,
                                         SLEQP_SOLVER_STATE_REAL_MAX_RAYLEIGH,
                                         &max_rayleigh));

  const double elapsed_seconds = sleqp_solver_get_elapsed_seconds(solver);

  double violation;

  SLEQP_CALL(sleqp_iterate_feasibility_residuum(problem,
                                                iterate,
                                                &violation));

  fprintf(stdout,
          "%s;%d;%d;%s;%f;%f;%d;%f;%d;%f;%f;%f\n",
          probname,
          num_variables,
          num_constraints,
          descriptions[status],
          sleqp_iterate_get_func_val(iterate),
          violation,
          iterations,
          elapsed_seconds,
          last_step_bdry,
          last_trust_radius,
          min_rayleigh,
          max_rayleigh);

  return SLEQP_OKAY;
}


int sleqp_cutest_run(const char* filename,
                     const char* probname,
                     const SleqpCutestOptions* cutest_options)
{
  integer funit = 42;        /* FORTRAN unit number for OUTSDIF.d */
  integer ierr;              /* Exit flag from OPEN and CLOSE */
  integer cutest_status;     /* Exit flag from CUTEst tools */


  integer CUTEst_nvar;        /* number of variables */
  integer CUTEst_ncons;       /* number of constraints */

  bool CUTest_constrained = false;

  ierr = 0;
  FORTRAN_open(&funit, filename, &ierr);

  if(ierr)
  {
    sleqp_log_error("Failed to open %s, aborting.", filename);
    return 1;
  }

  CUTEST_cdimen(&cutest_status, &funit, &CUTEst_nvar, &CUTEst_ncons);

  sleqp_log_info("Problem has %d variables, %d constraints",
                 CUTEst_nvar,
                 CUTEst_ncons);

  if(CUTEst_ncons)
  {
    CUTest_constrained = true;
  }

  SleqpCutestData* cutest_data;
  SleqpParams* params;

  SLEQP_CALL(sleqp_cutest_data_create(&cutest_data,
                                      funit,
                                      CUTEst_nvar,
                                      CUTEst_ncons));

  SLEQP_CALL(sleqp_params_create(&params));

  const double zero_eps = sleqp_params_get(params, SLEQP_PARAM_ZERO_EPS);

  SleqpSparseVec* x;

  SleqpOptions* options;
  SleqpProblem* problem;
  SleqpSolver* solver;

  SLEQP_CALL(sleqp_sparse_vector_create(&x, CUTEst_nvar, 0));

  SLEQP_CALL(sleqp_sparse_vector_from_raw(x, cutest_data->x, CUTEst_nvar, zero_eps));

  if(CUTest_constrained)
  {
    SLEQP_CALL(sleqp_cutest_cons_problem_create(&problem,
                                                cutest_data,
                                                params,
                                                cutest_options->force_nonlinear_constraints));
  }
  else
  {
    SLEQP_CALL(sleqp_cutest_uncons_problem_create(&problem,
                                                  cutest_data,
                                                  params));
  }

  SLEQP_CALL(sleqp_options_create(&options));

  if(cutest_options->enable_preprocessing)
  {
    SLEQP_CALL(sleqp_options_set_bool(options,
                                      SLEQP_OPTION_BOOL_ENABLE_PREPROCESSOR,
                                      true));
  }

  if(cutest_options->max_num_threads != SLEQP_NONE)
  {
    SLEQP_CALL(sleqp_options_set_int(options,
                                     SLEQP_OPTION_INT_NUM_THREADS,
                                     cutest_options->max_num_threads));
  }

  /*
  SLEQP_CALL(sleqp_options_set_int(options,
                                   SLEQP_OPTION_INT_DERIV_CHECK,
                                   SLEQP_DERIV_CHECK_FIRST));
  */

  SLEQP_CALL(sleqp_solver_create(&solver,
                                 problem,
                                 params,
                                 options,
                                 x,
                                 NULL));

  const int max_num_iterations = -1;
  const double time_limit = 3600;

  bool success = true;

  SLEQP_RETCODE retcode = sleqp_solver_solve(solver,
                                             max_num_iterations,
                                             time_limit);

  if(retcode == SLEQP_OKAY)
  {
    SLEQP_CALL(report_result(solver, problem, probname));

    SLEQP_STATUS status = sleqp_solver_get_status(solver);

    if(status == SLEQP_INVALID)
    {
      success = false;
    }
  }
  else
  {
    success = false;
  }

  /**/

  SLEQP_CALL(sleqp_solver_release(&solver));

  SLEQP_CALL(sleqp_options_release(&options));

  SLEQP_CALL(sleqp_problem_release(&problem));

  SLEQP_CALL(sleqp_sparse_vector_free(&x));

  SLEQP_CALL(sleqp_params_release(&params));

  SLEQP_CALL(sleqp_cutest_data_free(&cutest_data));

  FORTRAN_close(&funit, &ierr);

  if(ierr)
  {
    sleqp_log_error("Error closing %s on unit %d.\n",
                    filename,
                    funit);
  }

  return !(success);
}
