#include "solver.h"

#include "error.h"

SLEQP_RETCODE
sleqp_solver_add_callback(SleqpSolver* solver,
                          SLEQP_SOLVER_EVENT solver_event,
                          void* callback_func,
                          void* callback_data)
{
  if (solver_event < 0 || solver_event >= SLEQP_SOLVER_NUM_EVENTS)
  {
    sleqp_raise(SLEQP_ILLEGAL_ARGUMENT,
                "Invalid event number %d",
                solver_event);
  }

  SLEQP_CALL(sleqp_callback_handler_add(solver->callback_handlers[solver_event],
                                        callback_func,
                                        callback_data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_solver_remove_callback(SleqpSolver* solver,
                             SLEQP_SOLVER_EVENT solver_event,
                             void* callback_func,
                             void* callback_data)
{
  if (solver_event < 0 || solver_event >= SLEQP_SOLVER_NUM_EVENTS)
  {
    sleqp_raise(SLEQP_ILLEGAL_ARGUMENT,
                "Invalid event number %d",
                solver_event);
  }

  SLEQP_CALL(
    sleqp_callback_handler_remove(solver->callback_handlers[solver_event],
                                  callback_func,
                                  callback_data));

  return SLEQP_OKAY;
}
