#ifndef SLEQP_UNCONSTRAINED_CAUCHY_H
#define SLEQP_UNCONSTRAINED_CAUCHY_H

#include "cauchy.h"
#include "problem.h"

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_unconstrained_cauchy_create(SleqpCauchy** star,
                                  SleqpProblem* problem,
                                  SleqpSettings* settings);

#endif /* SLEQP_UNCONSTRAINED_CAUCHY_H */
