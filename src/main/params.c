#include "params.h"

#include <assert.h>

#include "log.h"
#include "mem.h"

struct SleqpParams
{
  int refcount;

  double values[SLEQP_NUM_PARAMS];
};

#define ZERO_EPS_DEFAULT 1e-20
#define EPS_DEFAULT 1e-10
#define OBJ_LOWER_DEFAULT -1e20

#define DERIV_PERTURBATION_DEFAULT 1e-8
#define DERIV_TOL_DEFAULT 1e-4

#define CAUCHY_TAU_DEFAULT 0.5
#define CAUCHY_ETA_DEFAULT 0.1

#define LINESEARCH_TAU_DEFAULT 0.5
#define LINESEARCH_ETA_DEFAULT 1e-4
#define LINESEARCH_CUTOFF_DEFAULT 1e-20

#define FEASIBILITY_TOL_DEFAULT 1e-6
#define SLACKNESS_TOL_DEFAULT 1e-6
#define STATIONARITY_TOL_DEFAULT 1e-6

#define ACCEPTED_REDUCTION_DEFAULT 1e-8

#define DEADPOINT_BOUND_DEFAULT 1e-12

SLEQP_RETCODE
sleqp_params_create(SleqpParams** star)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpParams* params = *star;

  *params = (SleqpParams){0};

  params->refcount = 1;

  params->values[SLEQP_PARAM_ZERO_EPS]           = ZERO_EPS_DEFAULT;
  params->values[SLEQP_PARAM_EPS]                = EPS_DEFAULT;
  params->values[SLEQP_PARAM_OBJ_LOWER]          = OBJ_LOWER_DEFAULT;
  params->values[SLEQP_PARAM_DERIV_PERTURBATION] = DERIV_PERTURBATION_DEFAULT;
  params->values[SLEQP_PARAM_DERIV_TOL]          = DERIV_TOL_DEFAULT;
  params->values[SLEQP_PARAM_CAUCHY_TAU]         = CAUCHY_TAU_DEFAULT;
  params->values[SLEQP_PARAM_CAUCHY_ETA]         = CAUCHY_ETA_DEFAULT;
  params->values[SLEQP_PARAM_LINESEARCH_TAU]     = LINESEARCH_TAU_DEFAULT;
  params->values[SLEQP_PARAM_LINESEARCH_ETA]     = LINESEARCH_ETA_DEFAULT;
  params->values[SLEQP_PARAM_LINESEARCH_CUTOFF]  = LINESEARCH_CUTOFF_DEFAULT;
  params->values[SLEQP_PARAM_FEASIBILITY_TOL]    = FEASIBILITY_TOL_DEFAULT;
  params->values[SLEQP_PARAM_SLACKNESS_TOL]      = SLACKNESS_TOL_DEFAULT;
  params->values[SLEQP_PARAM_STATIONARITY_TOL]   = STATIONARITY_TOL_DEFAULT;
  params->values[SLEQP_PARAM_ACCEPTED_REDUCTION] = ACCEPTED_REDUCTION_DEFAULT;
  params->values[SLEQP_PARAM_DEADPOINT_BOUND]    = DEADPOINT_BOUND_DEFAULT;

  return SLEQP_OKAY;
}

double
sleqp_params_value(const SleqpParams* params, SLEQP_PARAM param)
{
  assert(param >= 0);
  assert(param < SLEQP_NUM_PARAMS);

  return params->values[param];
}

SLEQP_RETCODE
sleqp_params_set_value(SleqpParams* params, SLEQP_PARAM param, double value)
{
  assert(param >= 0);
  assert(param < SLEQP_NUM_PARAMS);

  if (value <= 0.)
  {
    sleqp_log_error("Wrong value for parameter %d: %f", param, value);

    return SLEQP_ILLEGAL_ARGUMENT;
  }

  params->values[param] = value;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
params_free(SleqpParams** star)
{
  SleqpParams* params = *star;

  if (!params)
  {
    return SLEQP_OKAY;
  }

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_params_capture(SleqpParams* params)
{
  ++params->refcount;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_params_release(SleqpParams** star)
{
  SleqpParams* params = *star;

  if (!params)
  {
    return SLEQP_OKAY;
  }

  if (--params->refcount == 0)
  {
    SLEQP_CALL(params_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
