#ifndef SLEQP_DUAL_ESTIMATION_MIXED_H
#define SLEQP_DUAL_ESTIMATION_MIXED_H

#include "aug_jac/aug_jac.h"
#include "cauchy/cauchy.h"
#include "dual_estimation.h"

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_dual_estimation_mixed_create(SleqpDualEstimation** star,
                                   SleqpProblem* problem,
                                   SleqpCauchy* cauchy,
                                   SleqpAugJac* aug_jacobian);

#endif /* SLEQP_DUAL_ESTIMATION_MIXED_H */
