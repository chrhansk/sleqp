#include "sleqp_cutest_driver.h"

#include <assert.h>
#include <errno.h>
#include <string.h>

#include "cmp.h"
#include "iterate.h"
#include "log.h"
#include "mem.h"
#include "options.h"

#include "sleqp_cutest_constrained.h"
#include "sleqp_cutest_data.h"
#include "sleqp_cutest_types.h"
#include "sleqp_cutest_unconstrained.h"

static SLEQP_RETCODE
report_result(SleqpSolver* solver,
              SleqpProblem* problem,
              const char* probname,
              FILE* output)
{
  const char* descriptions[] = {
    [SLEQP_STATUS_UNKNOWN]         = "unknown",
    [SLEQP_STATUS_RUNNING]         = "running",
    [SLEQP_STATUS_OPTIMAL]         = "optimal",
    [SLEQP_STATUS_UNBOUNDED]       = "unbounded",
    [SLEQP_STATUS_ABORT_ITER]      = "abort_iter_limit",
    [SLEQP_STATUS_ABORT_TIME]      = "abort_time_limit",
    [SLEQP_STATUS_ABORT_DEADPOINT] = "abort_dead_point",
  };

  const int num_variables   = sleqp_problem_num_vars(problem);
  const int num_constraints = sleqp_problem_num_cons(problem);

  double* cache;

  SLEQP_CALL(
    sleqp_alloc_array(&cache, SLEQP_MAX(num_variables, num_constraints)));

  SLEQP_STATUS status = sleqp_solver_status(solver);
  SleqpIterate* iterate;

  SLEQP_CALL(sleqp_solver_solution(solver, &iterate));

  const int iterations = sleqp_solver_iterations(solver);

  int last_step_bdry;

  SLEQP_CALL(sleqp_solver_int_state(solver,
                                    SLEQP_SOLVER_STATE_INT_LAST_STEP_ON_BDRY,
                                    &last_step_bdry));

  double last_trust_radius;

  SLEQP_CALL(sleqp_solver_real_state(solver,
                                     SLEQP_SOLVER_STATE_REAL_TRUST_RADIUS,
                                     &last_trust_radius));

  double min_rayleigh;

  SLEQP_CALL(sleqp_solver_real_state(solver,
                                     SLEQP_SOLVER_STATE_REAL_MIN_RAYLEIGH,
                                     &min_rayleigh));

  double max_rayleigh;

  SLEQP_CALL(sleqp_solver_real_state(solver,
                                     SLEQP_SOLVER_STATE_REAL_MAX_RAYLEIGH,
                                     &max_rayleigh));

  const double elapsed_seconds = sleqp_solver_elapsed_seconds(solver);

  double feas_res;

  SLEQP_CALL(sleqp_iterate_feasibility_residuum(problem, iterate, &feas_res));

  double stat_res;

  SLEQP_CALL(
    sleqp_iterate_stationarity_residuum(problem, iterate, cache, &stat_res));

  double slack_res;

  SLEQP_CALL(sleqp_iterate_slackness_residuum(problem, iterate, &slack_res));

  fprintf(output,
          "%s;%d;%d;%s;%f;%.14e;%.14e;%.14e;%d;%f;%d;%f;%f;%f\n",
          probname,
          num_variables,
          num_constraints,
          descriptions[status],
          sleqp_iterate_obj_val(iterate),
          feas_res,
          slack_res,
          stat_res,
          iterations,
          elapsed_seconds,
          last_step_bdry,
          last_trust_radius,
          min_rayleigh,
          max_rayleigh);

  sleqp_free(&cache);

  return SLEQP_OKAY;
}

int
sleqp_cutest_run(const char* filename,
                 const char* probname,
                 const SleqpCutestOptions* cutest_options)
{
  integer funit = 42;    /* FORTRAN unit number for OUTSDIF.d */
  integer ierr;          /* Exit flag from OPEN and CLOSE */
  integer cutest_status; /* Exit flag from CUTEst tools */

  integer CUTEst_nvar;  /* number of variables */
  integer CUTEst_ncons; /* number of constraints */

  bool CUTest_constrained = false;

  FILE* output;

  if (cutest_options->output)
  {
    output = fopen(cutest_options->output, "w");

    if (!output)
    {
      sleqp_log_error("Failed to open %s: %s, aborting.",
                      cutest_options->output,
                      strerror(errno));

      return EXIT_FAILURE;
    }
  }
  else
  {
    output = stdout;
  }

  ierr = 0;
  FORTRAN_open(&funit, filename, &ierr);

  if (ierr)
  {
    sleqp_log_error("Failed to open %s, aborting.", filename);
    return 1;
  }

  if (!cutest_options->enable_logging)
  {
    sleqp_log_set_level(SLEQP_LOG_ERROR);
  }

  CUTEST_cdimen(&cutest_status, &funit, &CUTEst_nvar, &CUTEst_ncons);

  sleqp_log_info("Problem has %d variables, %d constraints",
                 CUTEst_nvar,
                 CUTEst_ncons);

  if (CUTEst_ncons)
  {
    CUTest_constrained = true;
  }

  SleqpCutestData* cutest_data;
  SleqpParams* params;

  SLEQP_CALL(
    sleqp_cutest_data_create(&cutest_data, funit, CUTEst_nvar, CUTEst_ncons));

  SLEQP_CALL(sleqp_params_create(&params));

  const double zero_eps = sleqp_params_value(params, SLEQP_PARAM_ZERO_EPS);

  SleqpSparseVec* x;

  SleqpOptions* options;
  SleqpProblem* problem;
  SleqpSolver* solver;

  SLEQP_CALL(sleqp_sparse_vector_create(&x, CUTEst_nvar, 0));

  // SLEQP_CALL(sleqp_sparse_vector_from_raw(x, sol, CUTEst_nvar, zero_eps));
  SLEQP_CALL(
    sleqp_sparse_vector_from_raw(x, cutest_data->x, CUTEst_nvar, zero_eps));

  if (CUTest_constrained)
  {
    SLEQP_CALL(sleqp_cutest_cons_problem_create(
      &problem,
      cutest_data,
      params,
      cutest_options->force_nonlinear_constraints));
  }
  else
  {
    SLEQP_CALL(
      sleqp_cutest_uncons_problem_create(&problem, cutest_data, params));
  }

  SLEQP_CALL(sleqp_options_create(&options));

  if (cutest_options->enable_preprocessing)
  {
    SLEQP_CALL(
      sleqp_options_set_bool_value(options,
                                   SLEQP_OPTION_BOOL_ENABLE_PREPROCESSOR,
                                   true));
  }

  if (cutest_options->max_num_threads != SLEQP_NONE)
  {
    SLEQP_CALL(sleqp_options_set_int_value(options,
                                           SLEQP_OPTION_INT_NUM_THREADS,
                                           cutest_options->max_num_threads));
  }

  /*
  SLEQP_CALL(sleqp_options_set_int_value(options,
                                         SLEQP_OPTION_INT_DERIV_CHECK,
                                         SLEQP_DERIV_CHECK_FIRST |
                                         SLEQP_DERIV_CHECK_SECOND_EXHAUSTIVE));
  */

  SLEQP_CALL(sleqp_solver_create(&solver, problem, params, options, x, NULL));

  const int max_num_iterations = -1;
  const double time_limit      = cutest_options->time_limit;

  bool success = true;

  SLEQP_RETCODE retcode
    = sleqp_solver_solve(solver, max_num_iterations, time_limit);

  if (retcode == SLEQP_OKAY)
  {
    SLEQP_CALL(report_result(solver, problem, probname, output));
  }
  else
  {
    sleqp_log_error("Failed to solve problem %s", probname);
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

  if (ierr)
  {
    sleqp_log_error("Error closing %s on unit %d.\n", filename, funit);
  }

  if (cutest_options->output)
  {
    fclose(output);
  }

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
  // return !(success);
}
