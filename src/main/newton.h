#ifndef SLEQP_NEWTON_H
#define SLEQP_NEWTON_H

/**
 * @file newton.h
 * @brief Definition of functions used for the computation of Newton (aka EQP) steps.
 **/

#include "options.h"
#include "params.h"
#include "problem.h"
#include "iterate.h"
#include "timer.h"
#include "working_step.h"

#include "aug_jac/aug_jac.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpNewtonData SleqpNewtonData;

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_newton_data_create(SleqpNewtonData** star,
                                         SleqpProblem* problem,
                                         SleqpWorkingStep* step,
                                         SleqpParams* params,
                                         SleqpOptions* options);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_newton_set_time_limit(SleqpNewtonData* data,
                                            double time_limit);

  SleqpTimer* sleqp_newton_get_timer(SleqpNewtonData* data);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_newton_set_iterate(SleqpNewtonData* data,
                                         SleqpIterate* iterate,
                                         SleqpAugJac* jacobian,
                                         double trust_radius,
                                         double penalty_parameter);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_newton_objective(SleqpNewtonData* data,
                                       const SleqpSparseVec* multipliers,
                                       const SleqpSparseVec* direction,
                                       double* objective);


  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_newton_add_violated_multipliers(SleqpNewtonData* data,
                                                      SleqpSparseVec* multipliers);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_newton_compute_step(SleqpNewtonData* data,
                                          SleqpSparseVec* multipliers,
                                          SleqpSparseVec* newton_step);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_newton_current_rayleigh(SleqpNewtonData* data,
                                              double* min_rayleigh,
                                              double* max_rayleigh);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_newton_data_capture(SleqpNewtonData* data);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_newton_data_release(SleqpNewtonData** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_NEWTON_H */
