#include "solver.h"

#include "feas.h"
#include "restoration.h"

static SLEQP_RETCODE
on_restoration_solver_accepted_iterate(SleqpProblemSolver* problem_solver,
                                       SleqpIterate* iterate,
                                       SleqpIterate* trial_iterate,
                                       void* callback_data)
{
  SleqpSolver* solver = (SleqpSolver*)callback_data;

  assert(problem_solver == solver->restoration_problem_solver);

  SleqpProblem* problem             = solver->problem;
  SleqpProblem* restoration_problem = solver->restoration_problem;
  SleqpFunc* restoration_func       = sleqp_problem_func(restoration_problem);

  const double feas_eps
    = sleqp_params_value(solver->params, SLEQP_PARAM_FEAS_TOL);

  SleqpVec* cons_val;

  SLEQP_CALL(sleqp_restoration_func_cons_val(restoration_func, &cons_val));

  double feas_res;

  SLEQP_CALL(sleqp_max_violation(problem, cons_val, &feas_res));

  if (feas_res <= feas_eps)
  {
    SLEQP_CALL(sleqp_problem_solver_abort(problem_solver));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
create_restoration_primal(SleqpSolver* solver)
{
  SleqpIterate* iterate
    = sleqp_problem_solver_get_iterate(solver->problem_solver);

  SLEQP_CALL(
    sleqp_restoration_problem_transform(solver->problem,
                                        sleqp_iterate_primal(iterate),
                                        sleqp_iterate_cons_val(iterate),
                                        solver->restoration_primal));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
create_primal(SleqpSolver* solver)
{
  SleqpIterate* restoration_iterate
    = sleqp_problem_solver_get_iterate(solver->restoration_problem_solver);

  SLEQP_CALL(
    sleqp_restoration_problem_restore(solver->problem,
                                      sleqp_iterate_primal(restoration_iterate),
                                      solver->primal));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
create_restoration_solver(SleqpSolver* solver)
{
  if (solver->restoration_problem_solver)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_restoration_problem_create(&solver->restoration_problem,
                                              solver->params,
                                              solver->problem));

  const int num_variables = sleqp_problem_num_vars(solver->restoration_problem);

  SLEQP_CALL(
    sleqp_vec_create_empty(&solver->restoration_primal, num_variables));

  SLEQP_CALL(create_restoration_primal(solver));

  SLEQP_CALL(sleqp_problem_solver_create(&solver->restoration_problem_solver,
                                         SLEQP_SOLVER_PHASE_RESTORATION,
                                         solver->restoration_problem,
                                         solver->params,
                                         solver->options,
                                         solver->restoration_primal));

  SLEQP_CALL(sleqp_problem_solver_add_callback(
    solver->restoration_problem_solver,
    SLEQP_PROBLEM_SOLVER_EVENT_ACCEPTED_ITERATE,
    on_restoration_solver_accepted_iterate,
    (void*)solver));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_solver_toggle_phase(SleqpSolver* solver)
{
  SLEQP_SOLVER_PHASE current_phase = solver->solver_phase;

  SLEQP_CALL(create_restoration_solver(solver));

  if (current_phase == SLEQP_SOLVER_PHASE_OPTIMIZATION)
  {
    sleqp_log_info("Switching to restoration phase");

    SLEQP_CALL(create_restoration_primal(solver));

    SLEQP_CALL(
      sleqp_problem_solver_set_primal(solver->restoration_problem_solver,
                                      solver->restoration_primal));

    solver->solver_phase = SLEQP_SOLVER_PHASE_RESTORATION;
  }
  else
  {
    sleqp_log_info("Switching to optimization phase");

    SLEQP_CALL(create_primal(solver));

    SLEQP_CALL(
      sleqp_problem_solver_set_primal(solver->problem_solver, solver->primal));

    solver->solver_phase = SLEQP_SOLVER_PHASE_OPTIMIZATION;
  }

  return SLEQP_OKAY;
}
