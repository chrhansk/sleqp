#ifndef SLEQP_STANDARD_CAUCHY_H
#define SLEQP_STANDARD_CAUCHY_H

/**
 * @file standard_cauchy.h
 * @brief Definition of LP-based Cauchy solver.
 **/

#include "cauchy.h"
#include "iterate.h"
#include "options.h"
#include "params.h"

#include "lp/lpi.h"

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_standard_cauchy_create(SleqpCauchy** star,
                             SleqpProblem* problem,
                             SleqpParams* params,
                             SleqpOptions* options);

#endif /* SLEQP_STANDARD_CAUCHY_H */
