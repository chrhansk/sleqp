#ifndef SLEQP_BOX_CONSTRAINED_CAUCHY_H
#define SLEQP_BOX_CONSTRAINED_CAUCHY_H

#include "cauchy.h"
#include "params.h"
#include "problem.h"

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_box_constrained_cauchy_create(SleqpCauchy** star,
                                    SleqpProblem* problem,
                                    SleqpParams* params);

#endif /* SLEQP_BOX_CONSTRAINED_CAUCHY_H */
