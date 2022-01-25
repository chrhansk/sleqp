#ifndef SLEQP_DUAL_ESTIMATION_H
#define SLEQP_DUAL_ESTIMATION_H

/**
 * @file dual_estimation.h
 * @brief Definition of functions for the estimation of dual variables.
 *
 * We follow the following convention:
 *
 * The Lagrangian is defined as
 * \f$ L(x, \lambda, \mu) = f(x) + \langle \lambda, c(x) \rangle +
 *     \langle x, \mu \rangle \f$.
 *
 * As a result, the signs of active dual variables are
 * non-negative for constraints / variables at their upper bounds and
 * non-positive for constraints / variables at their lower bounds.
 *
 **/

#include "dual_estimation_types.h"

typedef struct SleqpDualEstimation SleqpDualEstimation;

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_dual_estimation_create(SleqpDualEstimation** star,
                             SleqpDualEstimationCallbacks* callbacks,
                             void* estimation_data);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_estimate_duals(SleqpDualEstimation* estimation,
                     const SleqpIterate* iterate,
                     SleqpSparseVec* cons_dual,
                     SleqpSparseVec* vars_dual);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_dual_estimation_capture(SleqpDualEstimation* estimation);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_dual_estimation_release(SleqpDualEstimation** star);

#endif /* SLEQP_DUAL_ESTIMATION_H */
