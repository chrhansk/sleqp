#ifndef SLEQP_NEWTON_H
#define SLEQP_NEWTON_H

/**
 * @file sleqp_newton.h
 * @brief Definition of functions used for the computation of Newton (aka EQP) steps.
 **/

#include "sleqp_aug_jacobian.h"
#include "sleqp_options.h"
#include "sleqp_params.h"
#include "sleqp_problem.h"
#include "sleqp_iterate.h"
#include "sleqp_timer.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpNewtonData SleqpNewtonData;

  SLEQP_RETCODE sleqp_newton_data_create(SleqpNewtonData** star,
                                         SleqpProblem* problem,
                                         SleqpParams* params,
                                         SleqpOptions* options);

  SLEQP_RETCODE sleqp_newton_set_time_limit(SleqpNewtonData* data,
                                            double time_limit);

  SleqpTimer* sleqp_newton_get_timer(SleqpNewtonData* data);

  SLEQP_RETCODE sleqp_newton_set_iterate(SleqpNewtonData* data,
                                         SleqpIterate* iterate,
                                         SleqpAugJacobian* jacobian,
                                         double trust_radius,
                                         double penalty_parameter);

  SLEQP_RETCODE sleqp_newton_compute_multipliers(SleqpNewtonData* data,
                                                 SleqpSparseVec* multipliers);

  SLEQP_RETCODE sleqp_newton_compute_step(SleqpNewtonData* data,
                                          SleqpSparseVec* multipliers,
                                          SleqpSparseVec* newton_step);

  SLEQP_RETCODE sleqp_newton_data_capture(SleqpNewtonData* data);

  SLEQP_RETCODE sleqp_newton_data_release(SleqpNewtonData** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_NEWTON_H */
