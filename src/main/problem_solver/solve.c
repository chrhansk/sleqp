#include "problem_solver.h"

#include "cmp.h"
#include "error.h"
#include "feas.h"
#include "func.h"

static bool
exhausted_time_limit(SleqpProblemSolver* solver)
{
  return sleqp_timer_exhausted_time_limit(solver->elapsed_timer,
                                          solver->time_limit);
}

static SLEQP_RETCODE
print_warning(SleqpProblemSolver* solver)
{
  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;

  SleqpFunc* func = sleqp_problem_func(problem);

  SLEQP_FUNC_TYPE func_type = sleqp_func_get_type(func);

  const SLEQP_DERIV_CHECK deriv_check
    = sleqp_options_enum_value(solver->options, SLEQP_OPTION_ENUM_DERIV_CHECK);

  const int hessian_check_flags
    = (SLEQP_DERIV_CHECK_SECOND_EXHAUSTIVE | SLEQP_DERIV_CHECK_SECOND_SIMPLE);

  const bool hessian_check = (deriv_check & hessian_check_flags);

  if (hessian_check)
  {
    const bool inexact_hessian
      = sleqp_func_has_flags(func, SLEQP_FUNC_HESS_INEXACT);

    if (inexact_hessian)
    {
      sleqp_log_warn("Enabled second order derivative check while using a "
                     "quasi-Newton method");
    }
  }

  if (func_type == SLEQP_FUNC_TYPE_DYNAMIC)
  {
    const int check_first_flags = SLEQP_DERIV_CHECK_FIRST;

    const bool check_first = (deriv_check & check_first_flags);

    if (check_first)
    {
      sleqp_log_warn("Enabled first order derivative check while using a "
                     "dynamic function");
    }
  }

  {
    double total_violation;

    SLEQP_CALL(sleqp_total_violation(problem,
                                     sleqp_iterate_cons_val(iterate),
                                     &total_violation));

    const double obj_val = sleqp_iterate_obj_val(iterate);

    if (total_violation > 10. * SLEQP_ABS(obj_val))
    {
      sleqp_log_warn("Problem is badly scaled, constraint violation %g "
                     "significantly exceeds function value of %g",
                     total_violation,
                     obj_val);
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_problem_solver_solve(SleqpProblemSolver* solver,
                           int max_num_iterations,
                           double time_limit,
                           bool abort_on_local_infeasibility)
{
  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;

  const int num_variables   = sleqp_problem_num_vars(problem);
  const int num_constraints = sleqp_problem_num_cons(problem);

  solver->abort_on_local_infeasibility = abort_on_local_infeasibility;
  solver->status                       = SLEQP_PROBLEM_SOLVER_STATUS_RUNNING;

  {
    SleqpMat* cons_jac = sleqp_iterate_cons_jac(iterate);

    sleqp_log_info("Solving a problem with %d variables, %d constraints, %d "
                   "Jacobian nonzeros",
                   num_variables,
                   num_constraints,
                   sleqp_mat_nnz(cons_jac));
  }

  // Warnings
  SLEQP_CALL(print_warning(solver));

  solver->status = SLEQP_PROBLEM_SOLVER_STATUS_RUNNING;

  solver->time_limit = time_limit;
  solver->abort_next = false;

  solver->elapsed_iterations     = 0;
  solver->num_accepted_steps     = 0;
  solver->num_soc_accepted_steps = 0;
  solver->num_rejected_steps     = 0;
  solver->num_failed_eqp_steps   = 0;
  solver->last_step_type         = SLEQP_STEPTYPE_NONE;

  SLEQP_CALL(sleqp_timer_reset(solver->elapsed_timer));

  const double deadpoint_bound
    = sleqp_params_value(solver->params, SLEQP_PARAM_DEADPOINT_BOUND);

  SLEQP_CALL(sleqp_problem_solver_print_header(solver));

  // main solving loop
  while (true)
  {
    if (exhausted_time_limit(solver))
    {
      sleqp_log_info("Exhausted time limit, terminating");
      solver->status = SLEQP_PROBLEM_SOLVER_STATUS_ABORT_TIME;
      break;
    }

    if (max_num_iterations != SLEQP_NONE
        && solver->elapsed_iterations >= max_num_iterations)
    {
      sleqp_log_info("Reached iteration limit, terminating");
      solver->status = SLEQP_PROBLEM_SOLVER_STATUS_ABORT_ITER;
      break;
    }

    if (solver->abort_next)
    {
      sleqp_log_info("Abortion requested, terminating");
      solver->status = SLEQP_PROBLEM_SOLVER_STATUS_ABORT_MANUAL;
      break;
    }

    SLEQP_CALL(sleqp_timer_start(solver->elapsed_timer));

    SLEQP_RETCODE solver_status
      = sleqp_problem_solver_perform_iteration(solver);

    SLEQP_CALL(sleqp_timer_stop(solver->elapsed_timer));

    if (solver_status == SLEQP_ABORT_TIME || exhausted_time_limit(solver))
    {
      sleqp_log_info("Exhausted time limit, terminating");
      solver->status = SLEQP_PROBLEM_SOLVER_STATUS_ABORT_TIME;
      break;
    }

    SLEQP_CALL(solver_status);

    if (solver->lp_trust_radius <= deadpoint_bound
        || solver->trust_radius <= deadpoint_bound)
    {
      sleqp_log_warn("Reached dead point");
      solver->status = SLEQP_PROBLEM_SOLVER_STATUS_ABORT_DEADPOINT;
      break;
    }

    if (solver->status != SLEQP_PROBLEM_SOLVER_STATUS_RUNNING)
    {
      break;
    }
  }

  assert(solver->status != SLEQP_PROBLEM_SOLVER_STATUS_RUNNING);

  return SLEQP_OKAY;
}
