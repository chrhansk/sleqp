#include "solver.h"

#include "cmp.h"

double sleqp_solver_remaining_time(SleqpSolver* solver)
{
  double time_limit = solver->time_limit;

  if(time_limit != SLEQP_NONE)
  {
    double remaining_time = time_limit - sleqp_timer_get_ttl(solver->elapsed_timer);

    remaining_time = SLEQP_MAX(remaining_time, 0.);

    return remaining_time;
  }

  return SLEQP_NONE;
}
