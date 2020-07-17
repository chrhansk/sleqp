#ifndef SLEQP_OPTIONS_H
#define SLEQP_OPTIONS_H

/**
 * @file sleqp_options.h
 * @brief Definition of options.
 **/

#include "sleqp_types.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpOptions SleqpOptions;

  SLEQP_RETCODE sleqp_options_create(SleqpOptions** star);

  SLEQP_DERIV_CHECK sleqp_options_get_deriv_check(const SleqpOptions* options);

  SLEQP_HESSIAN_EVAL sleqp_options_get_hessian_eval(const SleqpOptions* options);

  SLEQP_DUAL_ESTIMATION_TYPE sleqp_options_get_dual_estimation_type(const SleqpOptions* options);

  int sleqp_options_get_quasi_newton_num_iterates(const SleqpOptions* options);

  SLEQP_RETCODE sleqp_options_set_deriv_check(SleqpOptions* options,
                                              SLEQP_DERIV_CHECK value);

  SLEQP_RETCODE sleqp_options_set_hessian_eval(SleqpOptions* options,
                                               SLEQP_HESSIAN_EVAL value);

  SLEQP_RETCODE sleqp_options_set_dual_estimation_type(SleqpOptions* options,
                                                       SLEQP_DUAL_ESTIMATION_TYPE dual_estimation_type);

  SLEQP_RETCODE sleqp_options_set_quasi_newton_num_iterates(SleqpOptions* options,
                                                            int size);

  SLEQP_RETCODE sleqp_options_free(SleqpOptions** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_OPTIONS_H */
