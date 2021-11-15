#include "solver.h"

#include "restoration.h"

static SLEQP_RETCODE
create_restoration_primal(SleqpSolver* solver)
{
  SleqpIterate* iterate
    = sleqp_problem_solver_get_iterate(solver->problem_solver);

  SLEQP_CALL(
    sleqp_restoration_problem_transform(solver->problem,
                                        sleqp_iterate_get_primal(iterate),
                                        sleqp_iterate_get_cons_val(iterate),
                                        solver->restoration_primal));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
create_primal(SleqpSolver* solver)
{
  SleqpIterate* restoration_iterate
    = sleqp_problem_solver_get_iterate(solver->restoration_problem_solver);

  SLEQP_CALL(sleqp_restoration_problem_restore(
    solver->problem,
    sleqp_iterate_get_primal(restoration_iterate),
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

  const int num_variables
    = sleqp_problem_num_variables(solver->restoration_problem);

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->restoration_primal,
                                              num_variables));

  SLEQP_CALL(create_restoration_primal(solver));

  SLEQP_CALL(sleqp_problem_solver_create(&solver->restoration_problem_solver,
                                         SLEQP_SOLVER_PHASE_RESTORATION,
                                         solver->restoration_problem,
                                         solver->params,
                                         solver->options,
                                         solver->restoration_primal));

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
