#ifndef SLEQP_DUAL_ESTIMATION_LSQ_H
#define SLEQP_DUAL_ESTIMATION_LSQ_H

#include "aug_jac/aug_jac.h"
#include "dual_estimation.h"

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_dual_estimation_lsq_create(SleqpDualEstimation** star,
                                 SleqpProblem* problem,
                                 SleqpAugJac* aug_jacobian);

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_estimate_duals_lsq(SleqpDualEstimation* estimation,
                         const SleqpIterate* iterate,
                         SleqpVec* cons_dual,
                         SleqpVec* vars_dual,
                         int* num_clipped_vars,
                         int* num_clipped_cons);

#endif /* SLEQP_DUAL_ESTIMATION_LSQ_H */
