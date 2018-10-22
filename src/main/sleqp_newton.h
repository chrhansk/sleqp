#ifndef SLEQP_NEWTON_H
#define SLEQP_NEWTON_H

#include "sleqp_aug_jacobian.h"
#include "sleqp_problem.h"
#include "sleqp_iterate.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpNewtonData SleqpNewtonData;

  SLEQP_RETCODE sleqp_newton_data_create(SleqpNewtonData** star,
                                         SleqpProblem* problem);

  SLEQP_RETCODE sleqp_newton_data_free(SleqpNewtonData** star);

  SLEQP_RETCODE sleqp_newton_compute_step(SleqpNewtonData* data,
                                          SleqpIterate* iterate,
                                          SleqpAugJacobian* jacobian,
                                          double trust_radius,
                                          double penalty_parameter,
                                          SleqpSparseVec* newton_step);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_NEWTON_H */
