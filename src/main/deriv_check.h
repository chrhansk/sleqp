#ifndef SLEQP_DERIV_CHECK_H
#define SLEQP_DERIV_CHECK_H

/**
 * @file deriv_check.h
 * @brief Definition of the derivative checker.
 **/

#include "func.h"
#include "iterate.h"
#include "params.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpDerivCheckData SleqpDerivCheckData;

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_deriv_checker_create(SleqpDerivCheckData** star,
                                           SleqpProblem* problem,
                                           SleqpParams* params);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_deriv_check_perform(SleqpDerivCheckData* data,
                                          SleqpIterate* iterate,
                                          SLEQP_DERIV_CHECK flags);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_deriv_checker_free(SleqpDerivCheckData** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_DERIV_CHECK_H */
