#ifndef SLEQP_SOC_H
#define SLEQP_SOC_H

/**
 * @file soc.h
 * @brief Definition of SOC (second-order correction) functions.
 **/

#include "aug_jacobian.h"
#include "iterate.h"
#include "problem.h"
#include "util.h"


#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpSOC SleqpSOC;

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_soc_data_create(SleqpSOC** star,
                                      SleqpProblem* problem,
                                      SleqpParams* params);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_soc_compute_direction(SleqpSOC* soc_data,
                                            SleqpAugJacobian* augmented_jacobian,
                                            const SleqpIterate* iterate,
                                            const SleqpIterate* trial_iterate,
                                            SleqpSparseVec* soc_direction);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_soc_compute_trial_point(SleqpSOC* soc_data,
                                              SleqpAugJacobian* augmented_jacobian,
                                              const SleqpIterate* iterate,
                                              const SleqpSparseVec* trial_step,
                                              const SleqpIterate* trial_iterate,
                                              SleqpSparseVec* soc_trial_point,
                                              double* soc_step_norm);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_soc_data_capture(SleqpSOC* soc_data);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_soc_data_release(SleqpSOC** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_SOC_H */
