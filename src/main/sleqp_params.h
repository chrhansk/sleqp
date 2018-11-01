#ifndef SLEQP_PARAMS_H
#define SLEQP_PARAMS_H

#include "sleqp_types.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpParams SleqpParams;

  SLEQP_RETCODE sleqp_params_create(SleqpParams** star);

  double sleqp_params_get_eps(SleqpParams* params);

  double sleqp_params_get_deriv_pertubation(SleqpParams* params);

  double sleqp_params_get_deriv_tolerance(SleqpParams* params);

  double sleqp_params_get_cauchy_tau(SleqpParams* params);
  double sleqp_params_get_cauchy_eta(SleqpParams* params);

  double sleqp_params_get_linesearch_tau(SleqpParams* params);
  double sleqp_params_get_linesearch_eta(SleqpParams* params);
  double sleqp_params_get_linesearch_cutoff(SleqpParams* params);

  double sleqp_params_get_optimality_tol(SleqpParams* params);

  double sleqp_params_get_accepted_reduction(SleqpParams* params);

  SLEQP_RETCODE sleqp_params_free(SleqpParams** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_PARAMS_H */
