#ifndef SLEQP_PUB_PARAMS_H
#define SLEQP_PUB_PARAMS_H

/**
 * @file pub_params.h
 * @brief Definition of numerical parameters.
 **/

#include "sleqp/export.h"
#include "sleqp/pub_types.h"

typedef enum
{
  SLEQP_PARAM_ZERO_EPS = 0,
  SLEQP_PARAM_EPS,
  SLEQP_PARAM_OBJ_LOWER,
  SLEQP_PARAM_DERIV_PERTURBATION,
  SLEQP_PARAM_DERIV_TOL,
  SLEQP_PARAM_CAUCHY_TAU,
  SLEQP_PARAM_CAUCHY_ETA,
  SLEQP_PARAM_LINESEARCH_TAU,
  SLEQP_PARAM_LINESEARCH_ETA,
  SLEQP_PARAM_LINESEARCH_CUTOFF,
  SLEQP_PARAM_FEAS_TOL,
  SLEQP_PARAM_SLACK_TOL,
  SLEQP_PARAM_STAT_TOL,
  SLEQP_PARAM_ACCEPTED_REDUCTION,
  SLEQP_PARAM_DEADPOINT_BOUND,
  SLEQP_NUM_PARAMS
} SLEQP_PARAM;

typedef struct SleqpParams SleqpParams;

SLEQP_EXPORT SLEQP_NODISCARD SLEQP_RETCODE
sleqp_params_create(SleqpParams** star);

SLEQP_EXPORT const char*
sleqp_params_name(SLEQP_PARAM param);

SLEQP_EXPORT const char*
sleqp_params_desc(SLEQP_PARAM param);

SLEQP_EXPORT double
sleqp_params_value(const SleqpParams* params, SLEQP_PARAM param);

SLEQP_EXPORT SLEQP_NODISCARD SLEQP_RETCODE
sleqp_params_set_value(SleqpParams* params, SLEQP_PARAM param, double value);

SLEQP_EXPORT SLEQP_NODISCARD SLEQP_RETCODE
sleqp_params_capture(SleqpParams* params);

SLEQP_EXPORT SLEQP_NODISCARD SLEQP_RETCODE
sleqp_params_release(SleqpParams** star);

#endif /* SLEQP_PUB_PARAMS_H */
