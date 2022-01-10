#ifndef SLEQP_PUB_OPTIONS_H
#define SLEQP_PUB_OPTIONS_H

/**
 * @file pub_options.h
 * @brief Definition of options.
 **/

#include "sleqp/export.h"
#include "sleqp/pub_types.h"

typedef enum
{
  SLEQP_OPTION_ENUM_DERIV_CHECK = 0,
  SLEQP_OPTION_ENUM_HESS_EVAL,
  SLEQP_OPTION_ENUM_DUAL_ESTIMATION_TYPE,
  SLEQP_OPTION_ENUM_FLOAT_WARNING_FLAGS,
  SLEQP_OPTION_ENUM_FLOAT_ERROR_FLAGS,
  SLEQP_OPTION_ENUM_BFGS_SIZING,
  SLEQP_OPTION_ENUM_TR_SOLVER,
  SLEQP_OPTION_ENUM_POLISHING_TYPE,
  SLEQP_OPTION_ENUM_STEP_RULE,
  SLEQP_OPTION_ENUM_LINESEARCH,
  SLEQP_OPTION_ENUM_PARAMETRIC_CAUCHY,
  SLEQP_OPTION_ENUM_INITIAL_TR_CHOICE,
  SLEQP_NUM_ENUM_OPTIONS
} SLEQP_OPTION_ENUM;

typedef enum
{
  SLEQP_OPTION_INT_NUM_QUASI_NEWTON_ITERATES = 0,
  SLEQP_OPTION_INT_MAX_NEWTON_ITERATIONS,
  SLEQP_OPTION_INT_NUM_THREADS,
  SLEQP_NUM_INT_OPTIONS
} SLEQP_OPTION_INT;

typedef enum
{
  SLEQP_OPTION_BOOL_PERFORM_NEWTON_STEP = 0,
  SLEQP_OPTION_BOOL_PERFORM_SOC,
  SLEQP_OPTION_BOOL_USE_QUADRATIC_MODEL,
  SLEQP_OPTION_BOOL_ALWAYS_WARM_START_LP,
  SLEQP_OPTION_BOOL_ENABLE_RESTORATION_PHASE,
  SLEQP_OPTION_BOOL_ENABLE_PREPROCESSOR,
  SLEQP_NUM_BOOL_OPTIONS
} SLEQP_OPTION_BOOL;

typedef struct SleqpOptions SleqpOptions;

SLEQP_EXPORT const char*
sleqp_options_enum_name(SLEQP_OPTION_ENUM options);

SLEQP_EXPORT const char*
sleqp_options_int_name(SLEQP_OPTION_INT options);

SLEQP_EXPORT const char*
sleqp_options_bool_name(SLEQP_OPTION_BOOL options);

SLEQP_EXPORT SLEQP_NODISCARD SLEQP_RETCODE
sleqp_options_create(SleqpOptions** star);

SLEQP_EXPORT
int
sleqp_options_enum_value(const SleqpOptions* options, SLEQP_OPTION_ENUM option);

SLEQP_EXPORT SLEQP_NODISCARD SLEQP_RETCODE
sleqp_options_set_enum_value(SleqpOptions* options,
                             SLEQP_OPTION_ENUM option,
                             int value);

SLEQP_EXPORT
int
sleqp_options_int_value(const SleqpOptions* options, SLEQP_OPTION_INT option);

SLEQP_EXPORT SLEQP_NODISCARD SLEQP_RETCODE
sleqp_options_set_int_value(SleqpOptions* options,
                            SLEQP_OPTION_INT option,
                            int value);

SLEQP_EXPORT
bool
sleqp_options_bool_value(const SleqpOptions* options, SLEQP_OPTION_BOOL option);

SLEQP_EXPORT SLEQP_NODISCARD SLEQP_RETCODE
sleqp_options_set_bool_value(SleqpOptions* options,
                             SLEQP_OPTION_BOOL option,
                             bool value);

SLEQP_EXPORT SLEQP_NODISCARD SLEQP_RETCODE
sleqp_options_capture(SleqpOptions* options);

SLEQP_EXPORT SLEQP_NODISCARD SLEQP_RETCODE
sleqp_options_release(SleqpOptions** star);

#endif /* SLEQP_PUB_OPTIONS_H */
