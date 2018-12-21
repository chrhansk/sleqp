#include "sleqp_params.h"

#include "sleqp_mem.h"

struct SleqpParams
{
  double eps;
  double zero_eps;

  double deriv_pertubation;
  double deriv_tolerance;

  double cauchy_tau;
  double cauchy_eta;

  double linesearch_tau;
  double linesearch_eta;
  double linesearch_cutoff;

  double optimality_tol;

  double accepted_reduction;
};

#define ZERO_EPS_DEFAULT 1e-16
#define EPS_DEFAULT 1e-10

#define DERIV_PERTUBATION_DEFAULT 1e-8
#define DERIV_TOLERANCE_DEFAULT 1e-4

#define CAUCHY_TAU_DEFAULT 0.5
#define CAUCHY_ETA_DEFAULT 0.1

#define LINESEARCH_TAU_DEFAULT 0.5
#define LINESEARCH_ETA_DEFAULT 1e-4
#define LINESEARCH_CUTOFF_DEFAULT 1e-20

#define OPTIMALITY_TOL_DEFAULT 1e-6

#define ACCEPTED_REDUCTION_DEFAULT 1e-8

SLEQP_RETCODE sleqp_params_create(SleqpParams** star)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpParams* params = *star;

  params->zero_eps = ZERO_EPS_DEFAULT;
  params->eps = EPS_DEFAULT;

  params->deriv_pertubation = DERIV_PERTUBATION_DEFAULT;
  params->deriv_tolerance = DERIV_TOLERANCE_DEFAULT;

  params->cauchy_tau = CAUCHY_TAU_DEFAULT;
  params->cauchy_eta = CAUCHY_ETA_DEFAULT;

  params->linesearch_tau = LINESEARCH_TAU_DEFAULT;
  params->linesearch_eta = LINESEARCH_ETA_DEFAULT;
  params->linesearch_cutoff = LINESEARCH_CUTOFF_DEFAULT;

  params->optimality_tol = OPTIMALITY_TOL_DEFAULT;

  params->accepted_reduction = ACCEPTED_REDUCTION_DEFAULT;

  return SLEQP_OKAY;
}

double sleqp_params_get_zero_eps(SleqpParams* params)
{
  return params->zero_eps;
}

double sleqp_params_get_eps(SleqpParams* params)
{
  return params->eps;
}

double sleqp_params_get_deriv_pertubation(SleqpParams* params)
{
  return params->deriv_pertubation;
}

double sleqp_params_get_deriv_tolerance(SleqpParams* params)
{
  return params->deriv_tolerance;
}

double sleqp_params_get_cauchy_tau(SleqpParams* params)
{
  return params->cauchy_tau;
}

double sleqp_params_get_cauchy_eta(SleqpParams* params)
{
  return params->cauchy_eta;
}

double sleqp_params_get_linesearch_tau(SleqpParams* params)
{
  return params->linesearch_tau;
}

double sleqp_params_get_linesearch_eta(SleqpParams* params)
{
  return params->linesearch_eta;
}

double sleqp_params_get_linesearch_cutoff(SleqpParams* params)
{
  return params->linesearch_cutoff;
}

double sleqp_params_get_optimality_tol(SleqpParams* params)
{
  return params->optimality_tol;
}

double sleqp_params_get_accepted_reduction(SleqpParams* params)
{
  return params->accepted_reduction;
}

SLEQP_RETCODE sleqp_params_free(SleqpParams** star)
{
  SleqpParams* params = *star;

  if(!params)
  {
    return SLEQP_OKAY;
  }

  sleqp_free(star);

  return SLEQP_OKAY;
}
