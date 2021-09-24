#ifndef SLEQP_BOX_CONSTRAINED_CAUCHY_H
#define SLEQP_BOX_CONSTRAINED_CAUCHY_H

#include "cauchy.h"
#include "problem.h"
#include "params.h"

#ifdef __cplusplus
extern "C" {
#endif

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_box_constrained_cauchy_create(SleqpCauchy** star,
                                                    SleqpProblem* problem,
                                                    SleqpParams* params);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_BOX_CONSTRAINED_CAUCHY_H */
