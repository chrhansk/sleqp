#include "dual_estimation.h"

#include "fail.h"
#include "mem.h"

struct SleqpDualEstimation
{
  int refcount;

  SleqpDualEstimationCallbacks callbacks;
  void* estimation_data;
};

SLEQP_RETCODE
sleqp_dual_estimation_create(SleqpDualEstimation** star,
                             SleqpDualEstimationCallbacks* callbacks,
                             void* estimation_data)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpDualEstimation* estimation = *star;

  *estimation = (SleqpDualEstimation){0};

  estimation->refcount = 1;

  estimation->callbacks       = *callbacks;
  estimation->estimation_data = estimation_data;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_estimate_duals(SleqpDualEstimation* estimation,
                     const SleqpIterate* iterate,
                     SleqpVec* cons_dual,
                     SleqpVec* vars_dual)
{
  SLEQP_CALL(estimation->callbacks.estimate_duals(iterate,
                                                  cons_dual,
                                                  vars_dual,
                                                  estimation->estimation_data));

#ifndef NDEBUG

  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);

  bool supports_cons_dual, supports_vars_dual;

  SLEQP_CALL(sleqp_working_set_supports_cons_dual(working_set,
                                                  cons_dual,
                                                  &supports_cons_dual));

  sleqp_assert_msg(supports_cons_dual,
                   "Working set does not support estimated constraint duals");

  SLEQP_CALL(sleqp_working_set_supports_vars_dual(working_set,
                                                  vars_dual,
                                                  &supports_vars_dual));

  sleqp_assert_msg(supports_vars_dual,
                   "Working set does not support estimated variable duals");

#endif

  return SLEQP_OKAY;
}

void*
sleqp_dual_estimation_data(SleqpDualEstimation* estimation)
{
  return estimation->estimation_data;
}

SLEQP_RETCODE
sleqp_dual_estimation_capture(SleqpDualEstimation* estimation)
{
  ++estimation->refcount;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
dual_estimation_free(SleqpDualEstimation** star)
{
  SleqpDualEstimation* estimation = *star;

  SLEQP_CALL(
    estimation->callbacks.estimation_free(estimation->estimation_data));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_dual_estimation_release(SleqpDualEstimation** star)
{
  SleqpDualEstimation* estimation = *star;

  if (!estimation)
  {
    return SLEQP_OKAY;
  }

  if (--estimation->refcount == 0)
  {
    SLEQP_CALL(dual_estimation_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
