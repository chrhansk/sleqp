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
