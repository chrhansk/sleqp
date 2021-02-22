#ifndef SLEQP_DERIV_CHECK_H
#define SLEQP_DERIV_CHECK_H

/**
 * @file sleqp_deriv_check.h
 * @brief Definition of the derivative checker.
 **/

#include "sleqp_func.h"
#include "sleqp_iterate.h"
#include "sleqp_params.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpDerivCheckData SleqpDerivCheckData;

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_deriv_checker_create(SleqpDerivCheckData** star,
                                           SleqpProblem* problem,
                                           SleqpParams* params);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_deriv_check_first_order(SleqpDerivCheckData* data,
                                              SleqpIterate* iterate);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_deriv_check_second_order_exhaustive(SleqpDerivCheckData* data,
                                                          SleqpIterate* iterate);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_deriv_check_second_order_simple(SleqpDerivCheckData* data,
                                                      SleqpIterate* iterate);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_deriv_checker_free(SleqpDerivCheckData** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_DERIV_CHECK_H */
