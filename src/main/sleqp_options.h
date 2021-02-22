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

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_options_create(SleqpOptions** star);

  SLEQP_EXPORT
  bool sleqp_options_get_perform_newton_step(const SleqpOptions* options);

  SLEQP_EXPORT
  bool sleqp_options_get_perform_soc(const SleqpOptions* options);

  SLEQP_EXPORT
  bool sleqp_options_get_use_quadratic_model(const SleqpOptions* options);

  SLEQP_EXPORT
  SLEQP_DERIV_CHECK sleqp_options_get_deriv_check(const SleqpOptions* options);

  SLEQP_EXPORT
  SLEQP_HESSIAN_EVAL sleqp_options_get_hessian_eval(const SleqpOptions* options);

  SLEQP_EXPORT
  SLEQP_DUAL_ESTIMATION_TYPE sleqp_options_get_dual_estimation_type(const SleqpOptions* options);

  SLEQP_EXPORT
  int sleqp_options_get_quasi_newton_num_iterates(const SleqpOptions* options);

  SLEQP_EXPORT
  int sleqp_options_get_max_newton_iterations(const SleqpOptions* options);

  SLEQP_EXPORT
  int sleqp_options_get_float_warning_flags(const SleqpOptions* options);

  SLEQP_EXPORT
  int sleqp_options_get_float_error_flags(const SleqpOptions* options);

  SLEQP_EXPORT
  SLEQP_BFGS_SIZING sleqp_options_get_bfgs_sizing(const SleqpOptions* options);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_options_set_perform_newton_step(SleqpOptions* options, bool value);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_options_set_perform_soc(SleqpOptions* options, bool value);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_options_set_use_quadratic_model(SleqpOptions* options, bool value);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_options_set_deriv_check(SleqpOptions* options,
                                              SLEQP_DERIV_CHECK value);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_options_set_hessian_eval(SleqpOptions* options,
                                               SLEQP_HESSIAN_EVAL value);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_options_set_dual_estimation_type(SleqpOptions* options,
                                                       SLEQP_DUAL_ESTIMATION_TYPE dual_estimation_type);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_options_set_quasi_newton_num_iterates(SleqpOptions* options,
                                                            int size);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_options_set_max_newton_iterations(SleqpOptions* options, int iterations);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_options_set_float_warning_flags(SleqpOptions* options, int flags);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_options_set_float_error_flags(SleqpOptions* options, int flags);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_options_set_bfgs_sizing(SleqpOptions* options, SLEQP_BFGS_SIZING sizing);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_options_capture(SleqpOptions* options);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_options_release(SleqpOptions** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_OPTIONS_H */
