#include "params.h"

#include "error.h"
#include "fail.h"
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
#define LINESEARCH_CUTOFF_DEFAULT 1e-6

#define FEASIBILITY_TOL_DEFAULT 1e-6
#define SLACKNESS_TOL_DEFAULT 1e-6
#define STATIONARITY_TOL_DEFAULT 1e-6

#define ACCEPTED_REDUCTION_DEFAULT 1e-8

#define DEADPOINT_BOUND_DEFAULT 1e-12

typedef struct
{
  const char* name;
  const char* desc;
} ParamInfo;

const ParamInfo param_infos[SLEQP_NUM_PARAMS] = {
  [SLEQP_PARAM_ZERO_EPS]           = {.name = "zero_eps",
                                      .desc = "Cut-off bound for "
                                              "sparse vectors / matrices"},
  [SLEQP_PARAM_EPS]                = {.name = "eps",
                                      .desc = "Basic accuracy"
                                              " for comparisons"},
  [SLEQP_PARAM_OBJ_LOWER]          = {.name = "obj_lower",
                                      .desc = "Lower bound on objective below which the "
                                              "problem is assumed to be unbounded"},
  [SLEQP_PARAM_DERIV_PERTURBATION] = {.name = "deriv_perturbation",
                                      .desc = "Perturbation used for "
                                              "the finite-differences "
                                              "used for checking "
                                              "derivatives"},
  [SLEQP_PARAM_DERIV_TOL]          = {.name = "deriv_tol",
                                      .desc = "Tolerance used "
                                              "for derivative checks"},
  [SLEQP_PARAM_CAUCHY_TAU]         = {.name = "cauchy_tau",
                                      .desc = "Reduction factor used"
                                              " during Cauchy line search"},
  [SLEQP_PARAM_CAUCHY_ETA]         = {.name = "cauchy_eta",
                                      .desc = "Efficiency parameter used "
                                              "during Cauchy line search"},
  [SLEQP_PARAM_LINESEARCH_TAU]     = {.name = "linesearch_tau",
                                      .desc = "Reduction factor used"
                                              " during trial point line search"},
  [SLEQP_PARAM_LINESEARCH_ETA]     = {.name = "linesearch_eta",
                                      .desc = "Efficiency parameter used "
                                              "during trial line search"},
  [SLEQP_PARAM_LINESEARCH_CUTOFF]  = {.name = "linesearch_cutoff",
                                      .desc = "Cutoff used during trial point "
                                              "line search below which a "
                                              "step size of 0 is chosen"},
  [SLEQP_PARAM_FEAS_TOL]           = {.name = "feas_tol",
                                      .desc = "Tolerance used to "
                                              "determine feasibility"},
  [SLEQP_PARAM_SLACK_TOL]          = {.name = "slack_tol",
                                      .desc = "Tolerance used to ensure "
                                              "complementary slackness"},
  [SLEQP_PARAM_STAT_TOL]           = {.name = "stat_tol",
                                      .desc = "Tolerance used to "
                                              "ensure stationarity"},
  [SLEQP_PARAM_ACCEPTED_REDUCTION] = {.name = "accepted_reduction",
                                      .desc = "Minimum model reduction "
                                              "to be accepted"},
  [SLEQP_PARAM_DEADPOINT_BOUND]    = {.name = "deadpoint_bound",
                                      .desc = "Lower bound on the normal "
                                              "and LP trust radius below which a "
                                              "stall is assumed"},
};

SLEQP_EXPORT const char*
sleqp_params_name(SLEQP_PARAM param)
{
  assert(param >= 0);
  assert(param < SLEQP_NUM_PARAMS);
  return param_infos[param].name;
}

SLEQP_EXPORT const char*
sleqp_params_desc(SLEQP_PARAM param)
{
  assert(param >= 0);
  assert(param < SLEQP_NUM_PARAMS);
  return param_infos[param].desc;
}

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
  params->values[SLEQP_PARAM_FEAS_TOL]           = FEASIBILITY_TOL_DEFAULT;
  params->values[SLEQP_PARAM_SLACK_TOL]          = SLACKNESS_TOL_DEFAULT;
  params->values[SLEQP_PARAM_STAT_TOL]           = STATIONARITY_TOL_DEFAULT;
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
    sleqp_raise(SLEQP_ILLEGAL_ARGUMENT,
                "Wrong value for parameter %d: %f",
                param,
                value);
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
