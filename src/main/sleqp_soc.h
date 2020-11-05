#ifndef SLEQP_SOC_H
#define SLEQP_SOC_H

/**
 * @file sleqp_soc.h
 * @brief Definition of SOC (second-order correction) functions.
 **/

#include "sleqp_aug_jacobian.h"
#include "sleqp_iterate.h"
#include "sleqp_problem.h"


#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpSOCData SleqpSOCData;

  SLEQP_RETCODE sleqp_soc_data_create(SleqpSOCData** star,
                                      SleqpProblem* problem,
                                      SleqpParams* params);

  SLEQP_RETCODE sleqp_soc_compute(SleqpSOCData* soc_data,
                                  SleqpAugJacobian* augmented_jacobian,
                                  SleqpIterate* iterate,
                                  SleqpIterate* trial_iterate,
                                  SleqpSparseVec* soc_direction);

  SLEQP_RETCODE sleqp_soc_data_capture(SleqpSOCData* soc_data);

  SLEQP_RETCODE sleqp_soc_data_release(SleqpSOCData** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_SOC_H */
