#ifndef SLEQP_DUAL_ESTIMATION_LSQ_H
#define SLEQP_DUAL_ESTIMATION_LSQ_H

#include "aug_jac/aug_jac.h"
#include "dual_estimation.h"

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_dual_estimation_lsq_create(SleqpDualEstimation** star,
                                 SleqpProblem* problem,
                                 SleqpAugJac* aug_jacobian);

#endif /* SLEQP_DUAL_ESTIMATION_LSQ_H */
