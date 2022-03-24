#include "dual_estimation_lp.h"

static SLEQP_RETCODE
estimate_duals_lp(const SleqpIterate* iterate,
                  SleqpVec* cons_dual,
                  SleqpVec* vars_dual,
                  void* estimation_data)
{
  SleqpCauchy* cauchy = (SleqpCauchy*)estimation_data;

  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);

  SLEQP_CALL(
    sleqp_cauchy_estimate_duals(cauchy, working_set, cons_dual, vars_dual));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
estimation_free(void* estimation_data)
{
  SleqpCauchy* cauchy = (SleqpCauchy*)estimation_data;

  SLEQP_CALL(sleqp_cauchy_release(&cauchy));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_dual_estimation_lp_create(SleqpDualEstimation** star, SleqpCauchy* cauchy)
{
  SleqpDualEstimationCallbacks callbacks = {
    .estimate_duals  = estimate_duals_lp,
    .estimation_free = estimation_free,
  };

  SLEQP_CALL(sleqp_cauchy_capture(cauchy));

  SLEQP_CALL(sleqp_dual_estimation_create(star, &callbacks, (void*)cauchy));

  return SLEQP_OKAY;
}
