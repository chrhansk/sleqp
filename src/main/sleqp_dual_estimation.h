#ifndef SLEQP_DUAL_ESTIMATION_H
#define SLEQP_DUAL_ESTIMATION_H

#include "sleqp.h"
#include "sleqp_active_set.h"
#include "sleqp_iterate.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpDualEstimationData SleqpDualEstimationData;

  SLEQP_RETCODE sleqp_dual_estimation_data_create(SleqpDualEstimationData** star,
                                                  SleqpProblem* problem);

  SLEQP_RETCODE sleqp_dual_estimation_compute(SleqpDualEstimationData* data,
                                              SleqpIterate* iterate,
                                              SleqpActiveSet* active_set);

  SLEQP_RETCODE sleqp_dual_estimation_data_free(SleqpDualEstimationData** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_DUAL_ESTIMATION_H */
