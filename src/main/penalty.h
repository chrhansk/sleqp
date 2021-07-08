#ifndef SLEQP_PENALTY_H
#define SLEQP_PENALTY_H

#include "cauchy.h"
#include "problem.h"

#ifdef __cplusplus
extern "C" {
#endif

  SLEQP_RETCODE sleqp_update_penalty(SleqpProblem* problem,
                                     SleqpIterate* iterate,
                                     SleqpCauchy* cauchy_data,
                                     double* penalty_parameter);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_PENALTY_H */
