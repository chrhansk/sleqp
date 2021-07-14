#include "sleqp_cutest_driver.h"

#include <assert.h>

#include "iterate.h"
#include "log.h"
#include "mem.h"

#include "sleqp_cutest_defs.h"
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
                     const char* probname)
{
  integer funit = 42;        /* FORTRAN unit number for OUTSDIF.d */
  integer iout = 6;          /* FORTRAN unit number for error output */
  integer ierr;              /* Exit flag from OPEN and CLOSE */
  integer io_buffer = 11;    /* FORTRAN unit internal input/output */
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

  double* x_dense, *var_ub_dense, *var_lb_dense;
  double* cons_ub_dense = NULL, *cons_lb_dense  = NULL;

  double *v = NULL;

  logical *equatn = NULL, *linear = NULL;
  int num_general = 0, num_linear = 0;

  // l_order = 2: general linear constraints should follow the
  //              general nonlinear ones
  integer e_order = 0, l_order = 2, v_order = 0;

  SLEQP_CALL(sleqp_alloc_array(&x_dense, CUTEst_nvar));
  SLEQP_CALL(sleqp_alloc_array(&var_lb_dense, CUTEst_nvar));
  SLEQP_CALL(sleqp_alloc_array(&var_ub_dense, CUTEst_nvar));

  if(CUTest_constrained)
  {
    SLEQP_CALL(sleqp_alloc_array(&cons_lb_dense, CUTEst_ncons));
    SLEQP_CALL(sleqp_alloc_array(&cons_ub_dense, CUTEst_ncons));

    SLEQP_CALL(sleqp_alloc_array(&equatn, CUTEst_ncons));
    SLEQP_CALL(sleqp_alloc_array(&linear, CUTEst_ncons));

    SLEQP_CALL(sleqp_alloc_array(&v, CUTEst_ncons));

    CUTEST_csetup(&cutest_status, &funit, &iout, &io_buffer,
                  &CUTEst_nvar, &CUTEst_ncons,
                  x_dense,
                  var_lb_dense,
                  var_ub_dense,
                  v,
                  cons_lb_dense,
                  cons_ub_dense,
                  equatn,
                  linear,
                  &e_order, &l_order, &v_order);

    {
      int i = 0;
      for(; i < CUTEst_ncons; ++i, ++num_general)
      {
        if(linear[i])
        {
          break;
        }
      }

      // ensure that l_order = 2 is satisfied
      for(; i < CUTEst_ncons; ++i)
      {
        assert(linear[i]);
      }
    }

    num_linear = CUTEst_ncons - num_general;
  }
  else
  {
    CUTEST_usetup(&cutest_status, &funit, &iout, &io_buffer,
                  &CUTEst_nvar, x_dense, var_lb_dense, var_ub_dense);
  }

  {
    const double inf = sleqp_infinity();

    for(int i = 0; i < CUTEst_nvar; i++)
    {
      if(var_lb_dense[i] == -CUTE_INF)
      {
        var_lb_dense[i] = -inf;
      }
      if(var_ub_dense[i] == CUTE_INF)
      {
        var_ub_dense[i] = inf;
      }
    }

    for (int i = 0; i < CUTEst_ncons; i++)
    {
      if(cons_lb_dense[i] == -CUTE_INF)
      {
        cons_lb_dense[i] = -inf;
      }
      if(cons_ub_dense[i] == CUTE_INF)
      {
        cons_ub_dense[i] = inf;
      }
    }
  }

  SleqpParams* params;

  SLEQP_CALL(sleqp_params_create(&params));

  const double eps = sleqp_params_get(params, SLEQP_PARAM_EPS);
  const double zero_eps = sleqp_params_get(params, SLEQP_PARAM_ZERO_EPS);

  SleqpSparseVec* var_lb;
  SleqpSparseVec* var_ub;
  SleqpSparseVec* x;

  SleqpSparseVec* cons_lb;
  SleqpSparseVec* cons_ub;
  SleqpFunc* func;

  SleqpSparseMatrix* linear_coeffs;
  SleqpSparseVec* sparse_cache;
  SleqpSparseVec* linear_lb;
  SleqpSparseVec* linear_ub;
  SleqpSparseVec* linear_offset;

  SleqpOptions* options;
  SleqpProblem* problem;
  SleqpSolver* solver;

  SLEQP_CALL(sleqp_sparse_vector_create(&var_lb, CUTEst_nvar, 0));
  SLEQP_CALL(sleqp_sparse_vector_create(&var_ub, CUTEst_nvar, 0));
  SLEQP_CALL(sleqp_sparse_vector_create(&x, CUTEst_nvar, 0));

  SLEQP_CALL(sleqp_sparse_vector_create(&cons_lb, num_general, 0));
  SLEQP_CALL(sleqp_sparse_vector_create(&cons_ub, num_general, 0));

  SLEQP_CALL(sleqp_sparse_vector_create(&sparse_cache, num_linear, 0));

  SLEQP_CALL(sleqp_sparse_vector_create(&linear_lb, num_linear, 0));
  SLEQP_CALL(sleqp_sparse_vector_create(&linear_ub, num_linear, 0));

  SLEQP_CALL(sleqp_sparse_vector_create(&linear_offset, num_linear, 0));

  SLEQP_CALL(sleqp_sparse_vector_from_raw(var_lb, var_lb_dense, CUTEst_nvar, zero_eps));
  SLEQP_CALL(sleqp_sparse_vector_from_raw(var_ub, var_ub_dense, CUTEst_nvar, zero_eps));

  SLEQP_CALL(sleqp_sparse_vector_from_raw(x, x_dense, CUTEst_nvar, zero_eps));

  SLEQP_CALL(sleqp_sparse_vector_from_raw(cons_lb, cons_lb_dense, num_general, zero_eps));
  SLEQP_CALL(sleqp_sparse_vector_from_raw(cons_ub, cons_ub_dense, num_general, zero_eps));

  SLEQP_CALL(sleqp_sparse_matrix_create(&linear_coeffs, num_linear, CUTEst_nvar, 0));

  if(CUTest_constrained)
  {
    SLEQP_CALL(sleqp_cutest_cons_func_create(&func,
                                             CUTEst_nvar,
                                             CUTEst_ncons,
                                             num_linear,
                                             params));

    SLEQP_CALL(sleqp_cutest_eval_linear(func, linear_coeffs));

    if(num_linear != 0)
    {
      SLEQP_CALL(sleqp_cutest_linear_offset(func,
                                            linear_offset));

      SLEQP_CALL(sleqp_sparse_vector_from_raw(sparse_cache,
                                              cons_lb_dense + num_general,
                                              num_linear,
                                              zero_eps));

      SLEQP_CALL(sleqp_sparse_vector_add_scaled(sparse_cache,
                                                linear_offset,
                                                1.,
                                                -1.,
                                                zero_eps,
                                                linear_lb));

      SLEQP_CALL(sleqp_sparse_vector_from_raw(sparse_cache,
                                              cons_ub_dense + num_general,
                                              num_linear,
                                              zero_eps));

      SLEQP_CALL(sleqp_sparse_vector_add_scaled(sparse_cache,
                                                linear_offset,
                                                1.,
                                                -1.,
                                                zero_eps,
                                                linear_ub));
    }
  }
  else
  {
    SLEQP_CALL(sleqp_cutest_uncons_func_create(&func, CUTEst_nvar, eps));
  }

  SLEQP_CALL(sleqp_problem_create(&problem,
                                  func,
                                  params,
                                  var_lb,
                                  var_ub,
                                  cons_lb,
                                  cons_ub,
                                  linear_coeffs,
                                  linear_lb,
                                  linear_ub));

  SLEQP_CALL(sleqp_options_create(&options));

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

  SLEQP_RETCODE status = sleqp_solver_solve(solver,
                                            max_num_iterations,
                                            time_limit);

  if(status == SLEQP_OKAY)
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

  SLEQP_CALL(sleqp_func_release(&func));

  SLEQP_CALL(sleqp_sparse_vector_free(&linear_offset));

  SLEQP_CALL(sleqp_sparse_vector_free(&linear_ub));
  SLEQP_CALL(sleqp_sparse_vector_free(&linear_lb));

  SLEQP_CALL(sleqp_sparse_matrix_release(&linear_coeffs));

  SLEQP_CALL(sleqp_sparse_vector_free(&sparse_cache));

  SLEQP_CALL(sleqp_sparse_vector_free(&cons_ub));
  SLEQP_CALL(sleqp_sparse_vector_free(&cons_lb));

  SLEQP_CALL(sleqp_sparse_vector_free(&x));
  SLEQP_CALL(sleqp_sparse_vector_free(&var_ub));
  SLEQP_CALL(sleqp_sparse_vector_free(&var_lb));

  SLEQP_CALL(sleqp_params_release(&params));

  sleqp_free(&v);
  sleqp_free(&linear);
  sleqp_free(&equatn);

  sleqp_free(&cons_ub_dense);
  sleqp_free(&cons_lb_dense);

  sleqp_free(&cons_ub_dense);
  sleqp_free(&cons_lb_dense);

  sleqp_free(&var_ub_dense);
  sleqp_free(&var_lb_dense);
  sleqp_free(&x_dense);


  FORTRAN_close(&funit, &ierr);

  if(ierr)
  {
    sleqp_log_error("Error closing %s on unit %d.\n",
                    filename,
                    funit);
  }

  return !(success);
}
