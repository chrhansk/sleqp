#include "problem_solver.h"

SLEQP_RETCODE
sleqp_problem_solver_add_callback(SleqpProblemSolver* solver,
                                  SLEQP_PROBLEM_SOLVER_EVENT solver_event,
                                  void* callback_func,
                                  void* callback_data)
{
  assert(solver_event >= 0);
  assert(solver_event < SLEQP_PROBLEM_SOLVER_NUM_EVENTS);

  SLEQP_CALL(sleqp_callback_handler_add(solver->callback_handlers[solver_event],
                                        callback_func,
                                        callback_data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_problem_solver_remove_callback(SleqpProblemSolver* solver,
                                     SLEQP_PROBLEM_SOLVER_EVENT solver_event,
                                     void* callback_func,
                                     void* callback_data)
{
  assert(solver_event >= 0);
  assert(solver_event < SLEQP_PROBLEM_SOLVER_NUM_EVENTS);

  SLEQP_CALL(
    sleqp_callback_handler_remove(solver->callback_handlers[solver_event],
                                  callback_func,
                                  callback_data));

  return SLEQP_OKAY;
}
