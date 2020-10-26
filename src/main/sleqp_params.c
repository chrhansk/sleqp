#include "sleqp_params.h"

#include "sleqp_mem.h"

struct SleqpParams
{
  double eps;
  double zero_eps;

  double deriv_perturbation;
  double deriv_tolerance;

  double cauchy_tau;
  double cauchy_eta;

  double linesearch_tau;
  double linesearch_eta;
  double linesearch_cutoff;

  double optimality_tol;

  double accepted_reduction;

  double deadpoint_bound;

  double newton_relative_tolerance;
};

#define ZERO_EPS_DEFAULT 1e-20
#define EPS_DEFAULT 1e-10

#define DERIV_PERTURBATION_DEFAULT 1e-8
#define DERIV_TOLERANCE_DEFAULT 1e-4

#define CAUCHY_TAU_DEFAULT 0.5
#define CAUCHY_ETA_DEFAULT 0.1

#define LINESEARCH_TAU_DEFAULT 0.5
#define LINESEARCH_ETA_DEFAULT 1e-4
#define LINESEARCH_CUTOFF_DEFAULT 1e-20

#define OPTIMALITY_TOL_DEFAULT 1e-6

#define ACCEPTED_REDUCTION_DEFAULT 1e-8

#define DEADPOINT_BOUND_DEFAULT 1e-10

#define NEWTON_RELATIVE_TOLERANCE_DEFAULT 1e-6

SLEQP_RETCODE sleqp_params_create(SleqpParams** star)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpParams* params = *star;

  params->zero_eps = ZERO_EPS_DEFAULT;
  params->eps = EPS_DEFAULT;

  params->deriv_perturbation = DERIV_PERTURBATION_DEFAULT;
  params->deriv_tolerance = DERIV_TOLERANCE_DEFAULT;

  params->cauchy_tau = CAUCHY_TAU_DEFAULT;
  params->cauchy_eta = CAUCHY_ETA_DEFAULT;

  params->linesearch_tau = LINESEARCH_TAU_DEFAULT;
  params->linesearch_eta = LINESEARCH_ETA_DEFAULT;
  params->linesearch_cutoff = LINESEARCH_CUTOFF_DEFAULT;

  params->optimality_tol = OPTIMALITY_TOL_DEFAULT;

  params->accepted_reduction = ACCEPTED_REDUCTION_DEFAULT;

  params->deadpoint_bound = DEADPOINT_BOUND_DEFAULT;

  params->newton_relative_tolerance = NEWTON_RELATIVE_TOLERANCE_DEFAULT;

  return SLEQP_OKAY;
}

double sleqp_params_get_zero_eps(const SleqpParams* params)
{
  return params->zero_eps;
}

double sleqp_params_get_eps(const SleqpParams* params)
{
  return params->eps;
}

double sleqp_params_get_deriv_perturbation(const SleqpParams* params)
{
  return params->deriv_perturbation;
}

double sleqp_params_get_deriv_tolerance(const SleqpParams* params)
{
  return params->deriv_tolerance;
}

double sleqp_params_get_cauchy_tau(const SleqpParams* params)
{
  return params->cauchy_tau;
}

double sleqp_params_get_cauchy_eta(const SleqpParams* params)
{
  return params->cauchy_eta;
}

double sleqp_params_get_linesearch_tau(const SleqpParams* params)
{
  return params->linesearch_tau;
}

double sleqp_params_get_linesearch_eta(const SleqpParams* params)
{
  return params->linesearch_eta;
}

double sleqp_params_get_linesearch_cutoff(const SleqpParams* params)
{
  return params->linesearch_cutoff;
}

double sleqp_params_get_optimality_tolerance(const SleqpParams* params)
{
  return params->optimality_tol;
}

double sleqp_params_get_accepted_reduction(const SleqpParams* params)
{
  return params->accepted_reduction;
}

double sleqp_params_get_deadpoint_bound(const SleqpParams* params)
{
  return params->deadpoint_bound;
}

double sleqp_params_get_newton_relative_tolerance(const SleqpParams* params)
{
  return params->newton_relative_tolerance;
}


SLEQP_RETCODE sleqp_params_set_zero_eps(SleqpParams* params, double value)
{
  params->zero_eps = value;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_params_set_eps(SleqpParams* params, double value)
{
  params->eps = value;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_params_set_deriv_perturbation(SleqpParams* params, double value)
{
  params->deriv_perturbation = value;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_params_set_deriv_tolerance(SleqpParams* params, double value)
{
  params->deriv_tolerance = value;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_params_set_cauchy_tau(SleqpParams* params, double value)
{
  params->cauchy_tau = value;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_params_set_cauchy_eta(SleqpParams* params, double value)
{
  params->cauchy_eta = value;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_params_set_linesearch_tau(SleqpParams* params, double value)
{
  params->linesearch_tau = value;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_params_set_linesearch_eta(SleqpParams* params, double value)
{
  params->linesearch_eta = value;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_params_set_linesearch_cutoff(SleqpParams* params, double value)
{
  params->linesearch_cutoff = value;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_params_set_optimality_tolerance(SleqpParams* params, double value)
{
  params->optimality_tol = value;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_params_set_accepted_reduction(SleqpParams* params, double value)
{
  params->accepted_reduction = value;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_params_set_deadpoint_bound(SleqpParams* params, double value)
{
  params->deadpoint_bound = value;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_params_set_newton_relative_tolerance(SleqpParams* params,
                                                         double value)
{
  params->newton_relative_tolerance = value;

  return SLEQP_OKAY;
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
