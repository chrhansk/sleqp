#ifndef DUAL_ESTIMATION_LP_H
#define DUAL_ESTIMATION_LP_H

#include "cauchy/cauchy.h"
#include "dual_estimation.h"

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_dual_estimation_lp_create(SleqpDualEstimation** star,
                                SleqpCauchy* cauchy);

#endif /* DUAL_ESTIMATION_LP_H */
