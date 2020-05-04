#ifndef SLEQP_NEWTON_H
#define SLEQP_NEWTON_H

/**
 * @file sleqp_newton.h
 * @brief Definition of functions used for the computation of Newton (aka EQP) steps.
 **/

#include "sleqp_aug_jacobian.h"
#include "sleqp_problem.h"
#include "sleqp_iterate.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpNewtonData SleqpNewtonData;

  SLEQP_RETCODE sleqp_newton_data_create(SleqpNewtonData** star,
                                         SleqpProblem* problem,
                                         SleqpParams* params);

  SLEQP_RETCODE sleqp_newton_set_time_limit(SleqpNewtonData* data,
                                            double time_limit);

  SleqpTimer* sleqp_newton_get_solve_timer(SleqpNewtonData* data);

  SLEQP_RETCODE sleqp_newton_compute_step(SleqpNewtonData* data,
                                          SleqpIterate* iterate,
                                          SleqpAugJacobian* jacobian,
                                          double trust_radius,
                                          double penalty_parameter,
                                          SleqpSparseVec* newton_step);

  SLEQP_RETCODE sleqp_newton_data_free(SleqpNewtonData** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_NEWTON_H */
