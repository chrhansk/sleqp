#ifndef SLEQP_UNCONSTRAINED_CAUCHY_H
#define SLEQP_UNCONSTRAINED_CAUCHY_H

#include "cauchy.h"
#include "problem.h"
#include "params.h"

#ifdef __cplusplus
extern "C" {
#endif

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_unconstrained_cauchy_create(SleqpCauchy** star,
                                                  SleqpProblem* problem,
                                                  SleqpParams* params);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_UNCONSTRAINED_CAUCHY_H */
