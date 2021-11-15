#ifndef SLEQP_STANDARD_CAUCHY_H
#define SLEQP_STANDARD_CAUCHY_H

/**
 * @file cauchy.h
 * @brief Definition of Cauchy step-related functions.
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
                             SleqpOptions* options,
                             SleqpLPi* lp_interface);

#endif /* SLEQP_STANDARD_CAUCHY_H */
