#ifndef SLEQP_BOX_CONSTRAINED_CAUCHY_H
#define SLEQP_BOX_CONSTRAINED_CAUCHY_H

#include "cauchy.h"
#include "problem.h"

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_box_constrained_cauchy_create(SleqpCauchy** star,
                                    SleqpProblem* problem,
                                    SleqpSettings* settings);

#endif /* SLEQP_BOX_CONSTRAINED_CAUCHY_H */
