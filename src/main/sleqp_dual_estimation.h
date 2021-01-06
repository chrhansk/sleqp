#ifndef SLEQP_DUAL_ESTIMATION_H
#define SLEQP_DUAL_ESTIMATION_H

/**
 * @file sleqp_dual_estimation.h
 * @brief Definition of functions for the estimation of dual variables.
 *
 * We follow the following convention:
 *
 * The Lagrangian is defined as
 * \f$ L(x, \lambda, \mu) = f(x) + \langle \lambda, c \rangle + \langle 1, \mu \rangle \f$.
 *
 * As a result, the signs of active dual variables are
 * non-negative for constraints / variables at their upper bounds and
 * non-positive for constraints / variables at their lower bounds.
 *
 *

 **/

#include "sleqp_working_set.h"
#include "sleqp_aug_jacobian.h"
#include "sleqp_iterate.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpDualEstimationData SleqpDualEstimationData;

  SLEQP_RETCODE sleqp_dual_estimation_data_create(SleqpDualEstimationData** star,
                                                  SleqpProblem* problem);

  /**
   * Computes the estimation of the dual variables for the given iterate
   * and stores the estimated valeus in the corresponding fields.
   *
   * @param[in]     estimation_data  The required estimation data
   * @param[in,out] iterate          The given iterate
   * @param[out]    residuum         The (optional) residuum of the estimation
   * @param[in]     aug_jacobian     The augmented Jacobian
   *
   **/
  SLEQP_RETCODE sleqp_dual_estimation_compute(SleqpDualEstimationData* estimation_data,
                                              SleqpIterate* iterate,
                                              SleqpSparseVec* residuum,
                                              SleqpAugJacobian* aug_jacobian);

  SLEQP_RETCODE sleqp_dual_estimation_data_free(SleqpDualEstimationData** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_DUAL_ESTIMATION_H */
