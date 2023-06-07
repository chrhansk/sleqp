#include "solver.h"

#include "util.h"

static SLEQP_RETCODE
check_feasibility(SleqpSolver* solver, bool* feasible)
{
  SleqpProblemSolver* problem_solver = solver->problem_solver;
  SleqpProblem* problem              = solver->problem;

  SleqpIterate* iterate = sleqp_problem_solver_iterate(problem_solver);

  bool reject;

  SLEQP_CALL(sleqp_set_and_evaluate(problem,
                                    iterate,
                                    SLEQP_VALUE_REASON_TRYING_ITERATE,
                                    &reject));

  *feasible = true;

  if (reject)
  {
    sleqp_log_debug("Function rejected restoration solution");

    *feasible = false;

    return SLEQP_OKAY;
  }

  double feasibility_residuum;

  SLEQP_CALL(sleqp_iterate_feasibility_residuum(problem,
                                                iterate,
                                                &feasibility_residuum));

  const double feas_eps
    = sleqp_settings_real_value(solver->settings, SLEQP_SETTINGS_REAL_FEAS_TOL);

  *feasible
    = sleqp_iterate_is_feasible(iterate, feasibility_residuum, feas_eps);

  return SLEQP_OKAY;
}

static SleqpProblemSolver*
current_problem_solver(SleqpSolver* solver)
{
  if (solver->solver_phase == SLEQP_SOLVER_PHASE_RESTORATION)
  {
    return solver->restoration_problem_solver;
  }

  return solver->problem_solver;
}

static SLEQP_RETCODE
set_initial_iterate(SleqpSolver* solver)
{
  assert(solver->solver_phase == SLEQP_SOLVER_PHASE_OPTIMIZATION);

  SleqpProblemSolver* problem_solver = solver->problem_solver;
  SleqpProblem* problem              = solver->problem;

  SLEQP_CALL(sleqp_problem_solver_set_cons_weights(problem_solver));

  SleqpIterate* iterate = sleqp_problem_solver_iterate(problem_solver);

  SleqpVec* primal = sleqp_iterate_primal(iterate);

  const double zero_eps
    = sleqp_settings_real_value(solver->settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  SLEQP_CALL(sleqp_vec_clip(solver->primal,
                            sleqp_problem_vars_lb(problem),
                            sleqp_problem_vars_ub(problem),
                            zero_eps,
                            primal));

  bool reject_initial;

  SLEQP_CALL(sleqp_set_and_evaluate(problem,
                                    iterate,
                                    SLEQP_VALUE_REASON_INIT,
                                    &reject_initial));

  if (reject_initial)
  {
    sleqp_raise(SLEQP_INTERNAL_ERROR, "Function rejected initial solution");
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
run_solving_loop(SleqpSolver* solver, int max_num_iterations, double time_limit)
{
  solver->iterations = 0;

  SLEQP_CALL(set_initial_iterate(solver));

  int remaining_iterations = max_num_iterations;

  bool enable_restoration
    = sleqp_settings_bool_value(solver->settings,
                               SLEQP_SETTINGS_BOOL_ENABLE_RESTORATION_PHASE);

  const bool unlimited_iterations = (max_num_iterations == SLEQP_NONE);

  while (true)
  {
    bool continue_loop = true;

    assert(remaining_iterations >= 0 || unlimited_iterations);

    if (sleqp_timer_exhausted_time_limit(solver->elapsed_timer, time_limit))
    {
      solver->status = SLEQP_STATUS_ABORT_TIME;
      break;
    }
    else if (remaining_iterations == 0)
    {
      solver->status = SLEQP_STATUS_ABORT_ITER;
      break;
    }

    double remaining_time
      = sleqp_timer_remaining_time(solver->elapsed_timer, time_limit);

    SleqpProblemSolver* problem_solver = current_problem_solver(solver);

    SLEQP_CALL(
      sleqp_problem_solver_set_iteration(problem_solver, solver->iterations));

    SLEQP_CALL(sleqp_problem_solver_solve(problem_solver,
                                          remaining_iterations,
                                          remaining_time,
                                          enable_restoration));

    SLEQP_PROBLEM_SOLVER_STATUS status
      = sleqp_problem_solver_status(problem_solver);

    assert(status != SLEQP_PROBLEM_SOLVER_STATUS_UNKNOWN);
    assert(status != SLEQP_PROBLEM_SOLVER_STATUS_RUNNING);

    solver->iterations
      += sleqp_problem_solver_elapsed_iterations(problem_solver);

    if (!unlimited_iterations)
    {
      remaining_iterations
        -= sleqp_problem_solver_elapsed_iterations(problem_solver);
    }

    // Propagate abort codes
    switch (status)
    {
    case SLEQP_PROBLEM_SOLVER_STATUS_ABORT_ITER:
      solver->status = SLEQP_STATUS_ABORT_ITER;
      continue_loop  = false;
      break;
    case SLEQP_PROBLEM_SOLVER_STATUS_ABORT_TIME:
      solver->status = SLEQP_STATUS_ABORT_TIME;
      continue_loop  = false;
      break;
    case SLEQP_PROBLEM_SOLVER_STATUS_ABORT_DEADPOINT:
      solver->status = SLEQP_STATUS_ABORT_DEADPOINT;
      continue_loop  = false;
      break;
    default:
      break;
    }

    if (!continue_loop)
    {
      break;
    }

    if (solver->solver_phase == SLEQP_SOLVER_PHASE_OPTIMIZATION)
    {
      switch (status)
      {
      case SLEQP_PROBLEM_SOLVER_STATUS_OPTIMAL:
        solver->status = SLEQP_STATUS_OPTIMAL;
        continue_loop  = false;
        break;
      case SLEQP_PROBLEM_SOLVER_STATUS_ABORT_MANUAL:
        solver->status = SLEQP_STATUS_ABORT_MANUAL;
        continue_loop  = false;
        break;
      case SLEQP_PROBLEM_SOLVER_STATUS_UNBOUNDED:
        solver->status = SLEQP_STATUS_UNBOUNDED;
        continue_loop  = false;
        break;
      case SLEQP_PROBLEM_SOLVER_STATUS_LOCALLY_INFEASIBLE:
        if (enable_restoration)
        {
          SLEQP_CALL(sleqp_solver_toggle_phase(solver));
        }
        break;
      default:
        // This should not happen
        assert(0);
      }
    }
    else
    {
      assert(solver->solver_phase == SLEQP_SOLVER_PHASE_RESTORATION);

      SLEQP_CALL(sleqp_solver_toggle_phase(solver));

      bool feasible;

      // manual abort must have been triggered by
      // "on_restoration_solver_accepted_iterate"
      // callback after achieving feasibility
      if (status == SLEQP_PROBLEM_SOLVER_STATUS_ABORT_MANUAL)
      {
        feasible = true;
#if SLEQP_DEBUG
        bool actually_feasible;
        SLEQP_CALL(check_feasibility(solver, &actually_feasible));
        assert(actually_feasible);
#endif
      }
      else
      {
        assert(status == SLEQP_PROBLEM_SOLVER_STATUS_OPTIMAL);
        SLEQP_CALL(check_feasibility(solver, &feasible));
      }

      if (!feasible)
      {
        sleqp_log_warn("Restoration failed");

        enable_restoration = false;
      }
    }

    if (!continue_loop)
    {
      break;
    }
  }

  assert(solver->status != SLEQP_STATUS_UNKNOWN);

  if (solver->solver_phase == SLEQP_SOLVER_PHASE_RESTORATION)
  {
    SLEQP_CALL(sleqp_solver_toggle_phase(solver));
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_solver_solve(SleqpSolver* solver,
                   int max_num_iterations,
                   double time_limit)
{
  if (solver->status == SLEQP_STATUS_INFEASIBLE)
  {
    sleqp_log_debug("Problem is infeasible, aborting");
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_timer_start(solver->elapsed_timer));

  SLEQP_CALL(run_solving_loop(solver, max_num_iterations, time_limit));

  SLEQP_CALL(sleqp_solver_restore_original_iterate(solver));

  double violation;

  SleqpIterate* iterate = sleqp_problem_solver_iterate(solver->problem_solver);

  SLEQP_CALL(
    sleqp_iterate_feasibility_residuum(solver->problem, iterate, &violation));

  SLEQP_POLISHING_TYPE polishing_type
    = sleqp_settings_enum_value(solver->settings,
                               SLEQP_SETTINGS_ENUM_POLISHING_TYPE);

  SLEQP_CALL(
    sleqp_polishing_polish(solver->polishing, iterate, polishing_type));

  SLEQP_CALLBACK_EVENT(solver->callback_handlers,
                       SLEQP_SOLVER_EVENT_FINISHED,
                       SLEQP_FINISHED,
                       solver,
                       solver->original_iterate);

  SLEQP_CALL(sleqp_timer_stop(solver->elapsed_timer));

  SLEQP_CALL(sleqp_solver_print_stats(solver, violation));

  return SLEQP_OKAY;
}
