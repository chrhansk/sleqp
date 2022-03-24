#include "solver.h"

SLEQP_RETCODE
sleqp_solver_real_state(const SleqpSolver* solver,
                        SLEQP_SOLVER_STATE_REAL state,
                        double* value)
{
  SLEQP_CALL(
    sleqp_problem_solver_get_real_state(solver->problem_solver, state, value));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_solver_int_state(const SleqpSolver* solver,
                       SLEQP_SOLVER_STATE_INT state,
                       int* value)
{
  SLEQP_CALL(
    sleqp_problem_solver_get_int_state(solver->problem_solver, state, value));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_solver_vec_state(const SleqpSolver* solver,
                       SLEQP_SOLVER_STATE_VEC value,
                       SleqpVec* result)
{
  SLEQP_CALL(
    sleqp_problem_solver_get_vec_state(solver->problem_solver, value, result));

  return SLEQP_OKAY;
}
