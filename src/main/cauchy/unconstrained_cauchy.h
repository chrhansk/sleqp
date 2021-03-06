#ifndef SLEQP_UNCONSTRAINED_CAUCHY_H
#define SLEQP_UNCONSTRAINED_CAUCHY_H

#include "cauchy.h"
#include "params.h"
#include "problem.h"

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_unconstrained_cauchy_create(SleqpCauchy** star,
                                  SleqpProblem* problem,
                                  SleqpParams* params);

#endif /* SLEQP_UNCONSTRAINED_CAUCHY_H */
