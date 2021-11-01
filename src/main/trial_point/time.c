#include "trial_point.h"

#include "cmp.h"


SLEQP_RETCODE
sleqp_trial_point_solver_set_time_limit(SleqpTrialPointSolver* solver,
                                        double time_limit)
{
  solver->time_limit = time_limit;

  SLEQP_CALL(sleqp_timer_reset(solver->elapsed_timer));

  return SLEQP_OKAY;
}

double sleqp_trial_point_solver_remaining_time(SleqpTrialPointSolver* solver)
{
  const double time_limit = solver->time_limit;

  if(time_limit != SLEQP_NONE)
  {
    double remaining_time = time_limit - sleqp_timer_get_ttl(solver->elapsed_timer);

    remaining_time = SLEQP_MAX(remaining_time, 0.);

    return remaining_time;
  }

  return SLEQP_NONE;
}
