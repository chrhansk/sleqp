#include "solver.h"

#include "log.h"
#include "lsq.h"

SLEQP_RETCODE
sleqp_solver_print_stats(SleqpSolver* solver, double violation)
{
  const char* descriptions[] = {
    [SLEQP_STATUS_UNKNOWN]
    = SLEQP_FORMAT_BOLD SLEQP_FORMAT_YELLOW "unknown" SLEQP_FORMAT_RESET,
    [SLEQP_STATUS_RUNNING]
    = SLEQP_FORMAT_BOLD SLEQP_FORMAT_YELLOW "running" SLEQP_FORMAT_RESET,
    [SLEQP_STATUS_OPTIMAL]
    = SLEQP_FORMAT_BOLD SLEQP_FORMAT_GREEN "optimal" SLEQP_FORMAT_RESET,
    [SLEQP_STATUS_INFEASIBLE]
    = SLEQP_FORMAT_BOLD SLEQP_FORMAT_RED "infeasible" SLEQP_FORMAT_RESET,
    [SLEQP_STATUS_UNBOUNDED]
    = SLEQP_FORMAT_BOLD SLEQP_FORMAT_RED "unbounded" SLEQP_FORMAT_RESET,
    [SLEQP_STATUS_ABORT_DEADPOINT] = SLEQP_FORMAT_BOLD SLEQP_FORMAT_RED
    "reached dead point" SLEQP_FORMAT_RESET,
    [SLEQP_STATUS_ABORT_ITER] = SLEQP_FORMAT_BOLD SLEQP_FORMAT_RED
    "reached iteration limit" SLEQP_FORMAT_RESET,
    [SLEQP_STATUS_ABORT_MANUAL]
    = SLEQP_FORMAT_BOLD SLEQP_FORMAT_RED "manual abort" SLEQP_FORMAT_RESET,
    [SLEQP_STATUS_ABORT_TIME] = SLEQP_FORMAT_BOLD SLEQP_FORMAT_RED
    "reached time limit" SLEQP_FORMAT_RESET,
  };

  SleqpIterate* iterate = sleqp_problem_solver_iterate(solver->problem_solver);

  SleqpFunc* original_func = sleqp_problem_func(solver->original_problem);
  SleqpFunc* func          = sleqp_problem_func(solver->problem);

  const bool with_hessian = !(solver->quasi_newton);

  double elapsed_seconds = sleqp_timer_get_ttl(solver->elapsed_timer);

  SleqpTimer* preprocessing_timer = NULL;

  if (solver->preprocessor)
  {
    preprocessing_timer = sleqp_preprocessor_get_timer(solver->preprocessor);

    double preprocessing_seconds = sleqp_timer_get_ttl(preprocessing_timer);

    elapsed_seconds += preprocessing_seconds;
  }

  sleqp_log_info(SLEQP_FORMAT_BOLD "%30s: %s" SLEQP_FORMAT_RESET,
                 "Solution status",
                 descriptions[solver->status]);

  if (solver->scaling_data)
  {
    double unscaled_violation;

    SLEQP_CALL(sleqp_iterate_feasibility_residuum(solver->original_problem,
                                                  solver->original_iterate,
                                                  &unscaled_violation));

    sleqp_log_info(SLEQP_FORMAT_BOLD "%30s:     %5.10e" SLEQP_FORMAT_RESET,
                   "Scaled objective value",
                   sleqp_iterate_obj_val(iterate));

    sleqp_log_info(SLEQP_FORMAT_BOLD "%30s:     %5.10e" SLEQP_FORMAT_RESET,
                   "Scaled violation",
                   violation);

    sleqp_log_info("%30s:     %5.10e",
                   "Original objective value",
                   sleqp_iterate_obj_val(solver->original_iterate));

    sleqp_log_info("%30s:     %5.10e",
                   "Original violation",
                   unscaled_violation);
  }
  else
  {
    sleqp_log_info(SLEQP_FORMAT_BOLD "%30s:     %5.10e" SLEQP_FORMAT_RESET,
                   "Objective value",
                   sleqp_iterate_obj_val(iterate));

    sleqp_log_info(SLEQP_FORMAT_BOLD "%30s:     %5.10e" SLEQP_FORMAT_RESET,
                   "Violation",
                   violation);
  }

  sleqp_log_info("%30s: %5d", "Iterations", sleqp_solver_iterations(solver));

  sleqp_log_info("%30s: %8.2fs", "Solving time", elapsed_seconds);

  if (preprocessing_timer)
  {
    SLEQP_CALL(sleqp_timer_display(preprocessing_timer,
                                   "Preprocessing",
                                   elapsed_seconds));
  }

  SLEQP_CALL(sleqp_timer_display(sleqp_func_get_set_timer(original_func),
                                 "Setting function values",
                                 elapsed_seconds));

  if (sleqp_func_get_type(original_func) == SLEQP_FUNC_TYPE_LSQ)
  {
    SLEQP_CALL(sleqp_timer_display(sleqp_lsq_func_residual_timer(original_func),
                                   "Residual evaluations",
                                   elapsed_seconds));

    SLEQP_CALL(sleqp_timer_display(sleqp_lsq_func_forward_timer(original_func),
                                   "Residual forward sweeps",
                                   elapsed_seconds));

    SLEQP_CALL(sleqp_timer_display(sleqp_lsq_func_adjoint_timer(original_func),
                                   "Residual adjoint sweeps",
                                   elapsed_seconds));
  }
  else
  {
    SLEQP_CALL(sleqp_timer_display(sleqp_func_get_val_timer(original_func),
                                   "Objective evaluations",
                                   elapsed_seconds));

    SLEQP_CALL(sleqp_timer_display(sleqp_func_get_grad_timer(original_func),
                                   "Gradient evaluations",
                                   elapsed_seconds));
  }

  SLEQP_CALL(sleqp_timer_display(sleqp_func_get_cons_val_timer(original_func),
                                 "Constraint evaluations",
                                 elapsed_seconds));

  SLEQP_CALL(sleqp_timer_display(sleqp_func_get_cons_jac_timer(original_func),
                                 "Jacobian evaluations",
                                 elapsed_seconds));

  if (with_hessian)
  {
    SLEQP_CALL(sleqp_timer_display(sleqp_func_get_hess_timer(original_func),
                                   "Hessian products",
                                   elapsed_seconds));
  }

  if (solver->quasi_newton)
  {
    SLEQP_CALL(sleqp_timer_display(sleqp_func_get_hess_timer(func),
                                   "quasi-Newton products",
                                   elapsed_seconds));

    SLEQP_CALL(
      sleqp_timer_display(sleqp_quasi_newton_update_timer(solver->quasi_newton),
                          "quasi-Newton updates",
                          elapsed_seconds));
  }

  if (solver->restoration_problem_solver)
  {
    sleqp_log_info(SLEQP_FORMAT_BOLD "%30s" SLEQP_FORMAT_RESET,
                   "Optimization phase");

    SLEQP_CALL(sleqp_problem_solver_print_stats(solver->problem_solver));

    sleqp_log_info(SLEQP_FORMAT_BOLD "%30s" SLEQP_FORMAT_RESET,
                   "Restoration phase");

    SLEQP_CALL(
      sleqp_problem_solver_print_stats(solver->restoration_problem_solver));
  }
  else
  {
    SLEQP_CALL(sleqp_problem_solver_print_stats(solver->problem_solver));
  }

  return SLEQP_OKAY;
}
