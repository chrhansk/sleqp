#ifndef SLEQP_PUB_PARAMS_H
#define SLEQP_PUB_PARAMS_H

/**
 * @file pub_params.h
 * @brief Definition of numerical parameters.
 **/

#include "sleqp/export.h"
#include "sleqp/pub_types.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef enum {
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
    SLEQP_PARAM_FEASIBILITY_TOL,
    SLEQP_PARAM_SLACKNESS_TOL,
    SLEQP_PARAM_STATIONARITY_TOL,
    SLEQP_PARAM_ACCEPTED_REDUCTION,
    SLEQP_PARAM_DEADPOINT_BOUND,
    SLEQP_NUM_PARAMS
  } SLEQP_PARAM;

  typedef struct SleqpParams SleqpParams;

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_params_create(SleqpParams** star);

  SLEQP_EXPORT double sleqp_params_get(const SleqpParams* params,
                                       SLEQP_PARAM param);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_params_set(SleqpParams* params,
                                 SLEQP_PARAM param,
                                 double value);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_params_capture(SleqpParams* params);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_params_release(SleqpParams** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_PUB_PARAMS_H */
