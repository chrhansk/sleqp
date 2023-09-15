#ifndef SLEQP_DERIV_CHECK_H
#define SLEQP_DERIV_CHECK_H

/**
 * @file deriv_check.h
 * @brief Definition of the derivative checker.
 **/

#include "func.h"
#include "iterate.h"

typedef struct SleqpDerivChecker SleqpDerivChecker;

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_deriv_checker_create(SleqpDerivChecker** star,
                           SleqpProblem* problem,
                           SleqpSettings* settings);

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_deriv_check_perform(SleqpDerivChecker* deriv_checker,
                          SleqpIterate* iterate,
                          SLEQP_DERIV_CHECK flags);

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_deriv_checker_free(SleqpDerivChecker** star);

#endif /* SLEQP_DERIV_CHECK_H */
