#ifndef SLEQP_PUB_OPTIONS_H
#define SLEQP_PUB_OPTIONS_H

/**
 * @file pub_options.h
 * @brief Definition of options.
 **/

#include "export.h"
#include "pub_types.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef enum {
    SLEQP_OPTION_INT_DERIV_CHECK = 0,
    SLEQP_OPTION_INT_HESSIAN_EVAL,
    SLEQP_OPTION_INT_DUAL_ESTIMATION_TYPE,
    SLEQP_OPTION_INT_NUM_QUASI_NEWTON_ITERATES,
    SLEQP_OPTION_INT_MAX_NEWTON_ITERATIONS,
    SLEQP_OPTION_INT_FLOAT_WARNING_FLAGS,
    SLEQP_OPTION_INT_FLOAT_ERROR_FLAGS,
    SLEQP_OPTION_INT_BFGS_SIZING,
    SLEQP_OPTION_INT_TR_SOLVER,
    SLEQP_OPTION_INT_LINESEARCH,
    SLEQP_OPTION_INT_PARAMETRIC_CAUCHY,
    SLEQP_NUM_INT_OPTIONS
  } SLEQP_OPTION_INT;

  typedef enum {
    SLEQP_OPTION_BOOL_PERFORM_NEWTON_STEP = 0,
    SLEQP_OPTION_BOOL_PERFORM_SOC,
    SLEQP_OPTION_BOOL_USE_QUADRATIC_MODEL,
    SLEQP_OPTION_BOOL_ALWAYS_WARM_START_LP,
    SLEQP_OPTION_BOOL_ENABLE_PREPROCESSOR,
    SLEQP_NUM_BOOL_OPTIONS
  } SLEQP_OPTION_BOOL;

  typedef struct SleqpOptions SleqpOptions;

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_options_create(SleqpOptions** star);

  SLEQP_EXPORT
  int sleqp_options_get_int(const SleqpOptions* options,
                            SLEQP_OPTION_INT option);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_options_set_int(SleqpOptions* options,
                                      SLEQP_OPTION_INT option,
                                      int value);

  SLEQP_EXPORT
  bool sleqp_options_get_bool(const SleqpOptions* options,
                              SLEQP_OPTION_BOOL option);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_options_set_bool(SleqpOptions* options,
                                       SLEQP_OPTION_BOOL option,
                                       bool value);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_options_capture(SleqpOptions* options);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_options_release(SleqpOptions** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_PUB_OPTIONS_H */
