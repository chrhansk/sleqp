#ifndef SLEQP_PARAMS_H
#define SLEQP_PARAMS_H

/**
 * @file sleqp_params.h
 * @brief Definition of numerical parameters.
 **/

#include "sleqp_types.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpParams SleqpParams;

  SLEQP_RETCODE sleqp_params_create(SleqpParams** star);

  double sleqp_params_get_zero_eps(SleqpParams* params);

  double sleqp_params_get_eps(SleqpParams* params);

  double sleqp_params_get_deriv_perturbation(SleqpParams* params);

  double sleqp_params_get_deriv_tolerance(SleqpParams* params);

  double sleqp_params_get_cauchy_tau(SleqpParams* params);
  double sleqp_params_get_cauchy_eta(SleqpParams* params);

  double sleqp_params_get_linesearch_tau(SleqpParams* params);
  double sleqp_params_get_linesearch_eta(SleqpParams* params);
  double sleqp_params_get_linesearch_cutoff(SleqpParams* params);

  double sleqp_params_get_optimality_tolerance(SleqpParams* params);

  double sleqp_params_get_accepted_reduction(SleqpParams* params);

  double sleqp_params_get_deadpoint_bound(SleqpParams* params);


  SLEQP_RETCODE sleqp_params_set_zero_eps(SleqpParams* params, double value);

  SLEQP_RETCODE sleqp_params_set_eps(SleqpParams* params, double value);

  SLEQP_RETCODE sleqp_params_set_deriv_perturbation(SleqpParams* params, double value);

  SLEQP_RETCODE sleqp_params_set_deriv_tolerance(SleqpParams* params, double value);

  SLEQP_RETCODE sleqp_params_set_cauchy_tau(SleqpParams* params, double value);
  SLEQP_RETCODE sleqp_params_set_cauchy_eta(SleqpParams* params, double value);

  SLEQP_RETCODE sleqp_params_set_linesearch_tau(SleqpParams* params, double value);
  SLEQP_RETCODE sleqp_params_set_linesearch_eta(SleqpParams* params, double value);
  SLEQP_RETCODE sleqp_params_set_linesearch_cutoff(SleqpParams* params, double value);

  SLEQP_RETCODE sleqp_params_set_optimality_tolerance(SleqpParams* params, double value);

  SLEQP_RETCODE sleqp_params_set_accepted_reduction(SleqpParams* params, double value);

  SLEQP_RETCODE sleqp_params_set_deadpoint_bound(SleqpParams* params, double value);


  SLEQP_RETCODE sleqp_params_free(SleqpParams** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_PARAMS_H */
