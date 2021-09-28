#ifndef SLEQP_DUAL_ESTIMATION_H
#define SLEQP_DUAL_ESTIMATION_H

/**
 * @file dual_estimation.h
 * @brief Definition of functions for the estimation of dual variables.
 *
 * We follow the following convention:
 *
 * The Lagrangian is defined as
 * \f$ L(x, \lambda, \mu) = f(x) + \langle \lambda, c \rangle + \langle x, \mu \rangle \f$.
 *
 * As a result, the signs of active dual variables are
 * non-negative for constraints / variables at their upper bounds and
 * non-positive for constraints / variables at their lower bounds.
 *
 *

 **/

#include "working_set.h"
#include "iterate.h"

#include "aug_jac/aug_jac.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpDualEstimation SleqpDualEstimation;

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_dual_estimation_create(SleqpDualEstimation** star,
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
  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_dual_estimation_compute(SleqpDualEstimation* estimation_data,
                                              SleqpIterate* iterate,
                                              SleqpSparseVec* residuum,
                                              SleqpAugJac* aug_jacobian);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_dual_estimation_free(SleqpDualEstimation** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_DUAL_ESTIMATION_H */
