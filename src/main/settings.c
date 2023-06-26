#include "settings.h"

#include <ctype.h>
#include <fenv.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "enum.h"
#include "error.h"
#include "fail.h"
#include "log.h"
#include "mem.h"
#include "pub_error.h"
#include "pub_log.h"
#include "pub_settings.h"
#include "pub_types.h"
#include "types.h"

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

#define PERFORM_NEWTON_DEFAULT true
#define GLOBAL_PENALTY_RESETS_DEFAULT true
#define REDUCED_AUG_JAC_DEFAULT false
#define PERFORM_SOC_DEFAULT true
#define USE_QUADRATIC_MODEL_DEFAULT true
#define ALWAYS_WARM_START_LP_DEFAULT true
#define ENABLE_PREPROCESSOR_DEFAULT false
#define ENABLE_RESTORATION_PHASE_DEFAULT true
#define LP_RESOLVES_DEFAULT true

#define DERIV_CHECK_DEFAULT SLEQP_DERIV_CHECK_SKIP
#define HESS_EVAL_DEFAULT SLEQP_HESS_EVAL_EXACT
#define DUAL_ESTIMATION_TYPE_DEFAULT SLEQP_DUAL_ESTIMATION_TYPE_LSQ
#define FLOAT_WARN_FLAGS_DEFAULT FE_ALL_EXCEPT
#define FLOAT_ERR_FLAGS_DEFAULT (FE_OVERFLOW | FE_DIVBYZERO | FE_INVALID)
#define BFGS_SIZING_DEFAULT SLEQP_BFGS_SIZING_CENTERED_OL
#define TR_SOLVER_DEFAULT SLEQP_TR_SOLVER_AUTO
#define POLISHING_TYPE_DEFAULT SLEQP_POLISHING_ZERO_DUAL
#define STEP_RULE_DEFAULT SLEQP_STEP_RULE_DIRECT
#define LINESEARCH_DEFAULT SLEQP_LINESEARCH_APPROX
#define PARAMETRIC_CAUCHY_DEFAULT SLEQP_PARAMETRIC_CAUCHY_DISABLED
#define INITIAL_TR_CHOICE_DEFAULT SLEQP_INITIAL_TR_CHOICE_NARROW
#define AUG_JAC_METHOD_DEFAULT SLEQP_AUG_JAC_AUTO

#define QUASI_NEWTON_SIZE_DEFAULT 5
#define MAX_NEWTON_ITERATIONS_DEFAULT 100
#define NUM_THREADS_DEFAULT SLEQP_NONE

#define CHECK_FLOAT_ENV                                                        \
  do                                                                           \
  {                                                                            \
    if (!(math_errhandling & MATH_ERREXCEPT))                                  \
    {                                                                          \
      sleqp_log_warn("Float point error handling is not supported, setting "   \
                     "options has no effect");                                 \
    }                                                                          \
  } while (false)

struct SleqpSettings
{
  int refcount;

  int enum_values[SLEQP_NUM_ENUM_SETTINGS];
  int int_values[SLEQP_NUM_INT_SETTINGS];
  bool bool_values[SLEQP_NUM_BOOL_SETTINGS];
  double real_values[SLEQP_NUM_REAL_SETTINGS];

  SLEQP_BFGS_SIZING bfgs_sizing;
  SLEQP_TR_SOLVER tr_solver;
};

typedef struct
{
  const char* name;
  const char* desc;
} OptionInfo;

OptionInfo enum_option_info[SLEQP_NUM_ENUM_SETTINGS] = {
  [SLEQP_SETTINGS_ENUM_DERIV_CHECK]
  = {.name = "deriv_check",
     .desc = "Which types of derivative check to perform"},
  [SLEQP_SETTINGS_ENUM_HESS_EVAL]
  = {.name = "hess_eval",
     .desc = "Whether or not to use a specific quasi-Newton method"},
  [SLEQP_SETTINGS_ENUM_DUAL_ESTIMATION_TYPE]
  = {.name = "dual_estimation_type", .desc = "How to estimate the dual values"},
  [SLEQP_SETTINGS_ENUM_FLOAT_WARNING_FLAGS]
  = {.name = "float_warning_flags",
     .desc = "Which floating point errors trigger warnings"},
  [SLEQP_SETTINGS_ENUM_FLOAT_ERROR_FLAGS]
  = {.name = "float_error_flags",
     .desc = "Which floating point errors trigger errors"},
  [SLEQP_SETTINGS_ENUM_BFGS_SIZING]
  = {.name = "bfgs_sizing", .desc = "How to size the BFGS method"},
  [SLEQP_SETTINGS_ENUM_TR_SOLVER]
  = {.name = "tr_solver", .desc = "Which trust-region solver to use"},
  [SLEQP_SETTINGS_ENUM_POLISHING_TYPE]
  = {.name = "polishing_type", .desc = "Which polishing type to use"},
  [SLEQP_SETTINGS_ENUM_STEP_RULE]
  = {.name = "step_rule", .desc = "Which step rule to use"},
  [SLEQP_SETTINGS_ENUM_LINESEARCH]
  = {.name = "linesearch",
     .desc = "Which line search to use for the trial point"},
  [SLEQP_SETTINGS_ENUM_PARAMETRIC_CAUCHY]
  = {.name = "parametric_cauchy",
     .desc = "Whether to use a parametric Cauchy line search"},
  [SLEQP_SETTINGS_ENUM_INITIAL_TR_CHOICE]
  = {.name = "initial_tr_choice",
     .desc = "How to chose the initial trust radius"},
  [SLEQP_SETTINGS_ENUM_AUG_JAC_METHOD]
  = {.name = "augmented Jacobian method",
     .desc = "How to solve the augmented Jacobian systems"},
};

const OptionInfo real_option_info[SLEQP_NUM_REAL_SETTINGS] = {
  [SLEQP_SETTINGS_REAL_ZERO_EPS] = {.name = "zero_eps",
                                    .desc = "Cut-off bound for "
                                            "sparse vectors / matrices"},
  [SLEQP_SETTINGS_REAL_EPS]      = {.name = "eps",
                                    .desc = "Basic accuracy"
                                                 " for comparisons"},
  [SLEQP_SETTINGS_REAL_OBJ_LOWER]
  = {.name = "obj_lower",
     .desc = "Lower bound on objective below which the "
             "problem is assumed to be unbounded"},
  [SLEQP_SETTINGS_REAL_DERIV_PERTURBATION] = {.name = "deriv_perturbation",
                                              .desc = "Perturbation used for "
                                                      "the finite-differences "
                                                      "used for checking "
                                                      "derivatives"},
  [SLEQP_SETTINGS_REAL_DERIV_TOL]          = {.name = "deriv_tol",
                                              .desc = "Tolerance used "
                                                               "for derivative checks"},
  [SLEQP_SETTINGS_REAL_CAUCHY_TAU]         = {.name = "cauchy_tau",
                                              .desc = "Reduction factor used"
                                                              " during Cauchy line search"},
  [SLEQP_SETTINGS_REAL_CAUCHY_ETA]         = {.name = "cauchy_eta",
                                              .desc = "Efficiency parameter used "
                                                              "during Cauchy line search"},
  [SLEQP_SETTINGS_REAL_LINESEARCH_TAU]
  = {.name = "linesearch_tau",
     .desc = "Reduction factor used"
             " during trial point line search"},
  [SLEQP_SETTINGS_REAL_LINESEARCH_ETA] = {.name = "linesearch_eta",
                                          .desc = "Efficiency parameter used "
                                                  "during trial line search"},
  [SLEQP_SETTINGS_REAL_LINESEARCH_CUTOFF]
  = {.name = "linesearch_cutoff",
     .desc = "Cutoff used during trial point "
             "line search below which a "
             "step size of 0 is chosen"},
  [SLEQP_SETTINGS_REAL_FEAS_TOL]           = {.name = "feas_tol",
                                              .desc = "Tolerance used to "
                                                                "determine feasibility"},
  [SLEQP_SETTINGS_REAL_SLACK_TOL]          = {.name = "slack_tol",
                                              .desc = "Tolerance used to ensure "
                                                               "complementary slackness"},
  [SLEQP_SETTINGS_REAL_STAT_TOL]           = {.name = "stat_tol",
                                              .desc = "Tolerance used to "
                                                                "ensure stationarity"},
  [SLEQP_SETTINGS_REAL_ACCEPTED_REDUCTION] = {.name = "accepted_reduction",
                                              .desc = "Minimum model reduction "
                                                      "to be accepted"},
  [SLEQP_SETTINGS_REAL_DEADPOINT_BOUND]
  = {.name = "deadpoint_bound",
     .desc = "Lower bound on the normal "
             "and LP trust radius below which a "
             "stall is assumed"},
};

const char*
sleqp_settings_real_name(SLEQP_SETTINGS_REAL param)
{
  assert(param >= 0);
  assert(param < SLEQP_NUM_REAL_SETTINGS);
  return real_option_info[param].name;
}

const char*
sleqp_settings_real_desc(SLEQP_SETTINGS_REAL param)
{
  assert(param >= 0);
  assert(param < SLEQP_NUM_REAL_SETTINGS);
  return real_option_info[param].desc;
}

const char*
sleqp_settings_enum_name(SLEQP_SETTINGS_ENUM option)
{
  assert(option >= 0);
  assert(option < SLEQP_NUM_ENUM_SETTINGS);
  return enum_option_info[option].name;
}

const char*
sleqp_settings_enum_desc(SLEQP_SETTINGS_ENUM option)
{
  assert(option >= 0);
  assert(option < SLEQP_NUM_ENUM_SETTINGS);
  return enum_option_info[option].desc;
}

OptionInfo int_option_info[SLEQP_NUM_INT_SETTINGS] = {
  [SLEQP_SETTINGS_INT_NUM_QUASI_NEWTON_ITERATES]
  = {.name = "num_quasi_newton_iterates",
     .desc = "Number of iterates to be stored"
             " for the quasi-Newton approximation "},
  [SLEQP_SETTINGS_INT_MAX_NEWTON_ITERATIONS]
  = {.name = "max_newton_iterations",
     .desc = "Maximum number of Newton iterations "
             "to be performed per iteration"},
  [SLEQP_SETTINGS_INT_NUM_THREADS]
  = {.name = "num_threads",
     .desc = "The maximum number of threads to be used."
             "Set to SLEQP_NONE to remove restriction."},
};

const char*
sleqp_settings_int_name(SLEQP_SETTINGS_INT option)
{
  assert(option >= 0);
  assert(option < SLEQP_NUM_INT_SETTINGS);
  return int_option_info[option].name;
}

const char*
sleqp_settings_int_desc(SLEQP_SETTINGS_INT option)
{
  assert(option >= 0);
  assert(option < SLEQP_NUM_INT_SETTINGS);
  return int_option_info[option].desc;
}

OptionInfo bool_option_info[SLEQP_NUM_BOOL_SETTINGS]
  = {[SLEQP_SETTINGS_BOOL_PERFORM_NEWTON_STEP]
     = {.name = "perform_newton_step",
        .desc = "Whether or not to perform Newton steps"},
     [SLEQP_SETTINGS_BOOL_GLOBAL_PENALTY_RESETS]
     = {.name = "global_penalty_resets",
        .desc = "Whether or not to perform global penalty resets"},
     [SLEQP_SETTINGS_BOOL_PERFORM_SOC]
     = {.name = "perform_soc",
        .desc = "Whether or not to perform a second-order correction"},
     [SLEQP_SETTINGS_BOOL_USE_QUADRATIC_MODEL]
     = {.name = "use_quadratic_model",
        .desc = "Whether to use a quadratic or linear model"},
     [SLEQP_SETTINGS_BOOL_ALWAYS_WARM_START_LP]
     = {.name = "always_warm_start_lp",
        .desc = "Whether to warm-start the LP from existing bases"},
     [SLEQP_SETTINGS_BOOL_ENABLE_RESTORATION_PHASE]
     = {.name = "enable_restoration_phase",
        .desc = "Whether to enable a restoration phase"},
     [SLEQP_SETTINGS_BOOL_ENABLE_PREPROCESSOR]
     = {.name = "enable_preprocessor",
        .desc = "Whether to enable the built-in preprocessor"},
     [SLEQP_SETTINGS_BOOL_LP_RESOLVES]
     = {.name = "lp_resolves",
        .desc = "Enable LP resolves in case of ambiguous optimal bases"}};

const char*
sleqp_settings_bool_name(SLEQP_SETTINGS_BOOL option)
{
  assert(option >= 0);
  assert(option < SLEQP_NUM_BOOL_SETTINGS);
  return bool_option_info[option].name;
}

const char*
sleqp_settings_bool_desc(SLEQP_SETTINGS_BOOL option)
{
  assert(option >= 0);
  assert(option < SLEQP_NUM_BOOL_SETTINGS);
  return bool_option_info[option].desc;
}

SLEQP_RETCODE
sleqp_settings_create(SleqpSettings** star)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpSettings* settings = *star;

  *settings = (SleqpSettings){0};

  *settings = (SleqpSettings){
    .refcount = 1,
    .enum_values
    = {[SLEQP_SETTINGS_ENUM_DERIV_CHECK] = DERIV_CHECK_DEFAULT,
       [SLEQP_SETTINGS_ENUM_HESS_EVAL]   = HESS_EVAL_DEFAULT,
       [SLEQP_SETTINGS_ENUM_DUAL_ESTIMATION_TYPE]
       = DUAL_ESTIMATION_TYPE_DEFAULT,
       [SLEQP_SETTINGS_ENUM_FLOAT_WARNING_FLAGS] = FLOAT_WARN_FLAGS_DEFAULT,
       [SLEQP_SETTINGS_ENUM_FLOAT_ERROR_FLAGS]   = FLOAT_ERR_FLAGS_DEFAULT,
       [SLEQP_SETTINGS_ENUM_BFGS_SIZING]         = BFGS_SIZING_DEFAULT,
       [SLEQP_SETTINGS_ENUM_TR_SOLVER]           = TR_SOLVER_DEFAULT,
       [SLEQP_SETTINGS_ENUM_POLISHING_TYPE]      = POLISHING_TYPE_DEFAULT,
       [SLEQP_SETTINGS_ENUM_STEP_RULE]           = STEP_RULE_DEFAULT,
       [SLEQP_SETTINGS_ENUM_LINESEARCH]          = LINESEARCH_DEFAULT,
       [SLEQP_SETTINGS_ENUM_PARAMETRIC_CAUCHY]   = PARAMETRIC_CAUCHY_DEFAULT,
       [SLEQP_SETTINGS_ENUM_INITIAL_TR_CHOICE]   = INITIAL_TR_CHOICE_DEFAULT,
       [SLEQP_SETTINGS_ENUM_AUG_JAC_METHOD]      = AUG_JAC_METHOD_DEFAULT},
    .int_values = {[SLEQP_SETTINGS_INT_NUM_QUASI_NEWTON_ITERATES]
                   = QUASI_NEWTON_SIZE_DEFAULT,
                   [SLEQP_SETTINGS_INT_MAX_NEWTON_ITERATIONS]
                   = MAX_NEWTON_ITERATIONS_DEFAULT,
                   [SLEQP_SETTINGS_INT_NUM_THREADS] = NUM_THREADS_DEFAULT},
    .bool_values
    = {[SLEQP_SETTINGS_BOOL_PERFORM_NEWTON_STEP] = PERFORM_NEWTON_DEFAULT,
       [SLEQP_SETTINGS_BOOL_GLOBAL_PENALTY_RESETS]
       = GLOBAL_PENALTY_RESETS_DEFAULT,
       [SLEQP_SETTINGS_BOOL_PERFORM_SOC]         = PERFORM_SOC_DEFAULT,
       [SLEQP_SETTINGS_BOOL_USE_QUADRATIC_MODEL] = USE_QUADRATIC_MODEL_DEFAULT,
       [SLEQP_SETTINGS_BOOL_ALWAYS_WARM_START_LP]
       = ALWAYS_WARM_START_LP_DEFAULT,
       [SLEQP_SETTINGS_BOOL_ENABLE_RESTORATION_PHASE]
       = ENABLE_RESTORATION_PHASE_DEFAULT,
       [SLEQP_SETTINGS_BOOL_ENABLE_PREPROCESSOR] = ENABLE_PREPROCESSOR_DEFAULT,
       [SLEQP_SETTINGS_BOOL_LP_RESOLVES]         = LP_RESOLVES_DEFAULT},
    .real_values
    = {[SLEQP_SETTINGS_REAL_ZERO_EPS]           = ZERO_EPS_DEFAULT,
       [SLEQP_SETTINGS_REAL_EPS]                = EPS_DEFAULT,
       [SLEQP_SETTINGS_REAL_OBJ_LOWER]          = OBJ_LOWER_DEFAULT,
       [SLEQP_SETTINGS_REAL_DERIV_PERTURBATION] = DERIV_PERTURBATION_DEFAULT,
       [SLEQP_SETTINGS_REAL_DERIV_TOL]          = DERIV_TOL_DEFAULT,
       [SLEQP_SETTINGS_REAL_CAUCHY_TAU]         = CAUCHY_TAU_DEFAULT,
       [SLEQP_SETTINGS_REAL_CAUCHY_ETA]         = CAUCHY_ETA_DEFAULT,
       [SLEQP_SETTINGS_REAL_LINESEARCH_TAU]     = LINESEARCH_TAU_DEFAULT,
       [SLEQP_SETTINGS_REAL_LINESEARCH_ETA]     = LINESEARCH_ETA_DEFAULT,
       [SLEQP_SETTINGS_REAL_LINESEARCH_CUTOFF]  = LINESEARCH_CUTOFF_DEFAULT,
       [SLEQP_SETTINGS_REAL_FEAS_TOL]           = FEASIBILITY_TOL_DEFAULT,
       [SLEQP_SETTINGS_REAL_SLACK_TOL]          = SLACKNESS_TOL_DEFAULT,
       [SLEQP_SETTINGS_REAL_STAT_TOL]           = STATIONARITY_TOL_DEFAULT,
       [SLEQP_SETTINGS_REAL_ACCEPTED_REDUCTION] = ACCEPTED_REDUCTION_DEFAULT,
       [SLEQP_SETTINGS_REAL_DEADPOINT_BOUND]    = DEADPOINT_BOUND_DEFAULT}};

  return SLEQP_OKAY;
}

static const SleqpEnum*
get_enum(SLEQP_SETTINGS_ENUM option)
{
  switch (option)
  {
  case SLEQP_SETTINGS_ENUM_DERIV_CHECK:
    return sleqp_enum_deriv_check();
  case SLEQP_SETTINGS_ENUM_HESS_EVAL:
    return sleqp_enum_hess_eval();
  case SLEQP_SETTINGS_ENUM_DUAL_ESTIMATION_TYPE:
    return sleqp_enum_dual_estimation();
  case SLEQP_SETTINGS_ENUM_FLOAT_WARNING_FLAGS:
  case SLEQP_SETTINGS_ENUM_FLOAT_ERROR_FLAGS:
    return NULL;
  case SLEQP_SETTINGS_ENUM_BFGS_SIZING:
    return sleqp_enum_bfgs_sizing();
  case SLEQP_SETTINGS_ENUM_TR_SOLVER:
    return sleqp_enum_tr_solver();
  case SLEQP_SETTINGS_ENUM_POLISHING_TYPE:
    return sleqp_enum_polishing_type();
  case SLEQP_SETTINGS_ENUM_STEP_RULE:
    return sleqp_enum_step_rule();
  case SLEQP_SETTINGS_ENUM_LINESEARCH:
    return sleqp_enum_linesearch();
  case SLEQP_SETTINGS_ENUM_PARAMETRIC_CAUCHY:
    return sleqp_enum_parametric_cauchy();
  case SLEQP_SETTINGS_ENUM_INITIAL_TR_CHOICE:
    return sleqp_enum_initial_tr();
  case SLEQP_SETTINGS_ENUM_AUG_JAC_METHOD:
    return sleqp_enum_aug_jac_method();
  default:
    assert(0);
  }
  return NULL;
}

static bool
valid_member(SLEQP_SETTINGS_ENUM option, int value)
{
  const SleqpEnum* settings_enum = get_enum(option);

  if (!settings_enum)
  {
    return true;
  }

  return sleqp_enum_member(settings_enum, value);
}

#define OPT_VALID(option, NUM_OPTIONS)                                         \
  (((option) >= 0) && ((option) < (NUM_OPTIONS)))

double
sleqp_settings_real_value(const SleqpSettings* settings,
                          SLEQP_SETTINGS_REAL param)
{
  if (!OPT_VALID(param, SLEQP_NUM_REAL_SETTINGS))
  {
    sleqp_raise(SLEQP_ILLEGAL_ARGUMENT, "Invalid real option (%d)", param);
  }

  return settings->real_values[param];
}

SLEQP_RETCODE
sleqp_settings_set_real_value(SleqpSettings* settings,
                              SLEQP_SETTINGS_REAL param,
                              double value)
{
  if (!OPT_VALID(param, SLEQP_NUM_REAL_SETTINGS))
  {
    sleqp_raise(SLEQP_ILLEGAL_ARGUMENT, "Invalid int option (%d)", param);
  }

  settings->real_values[param] = value;

  return SLEQP_OKAY;
}

int
sleqp_settings_enum_value(const SleqpSettings* settings,
                          SLEQP_SETTINGS_ENUM option)
{
  if (!OPT_VALID(option, SLEQP_NUM_ENUM_SETTINGS))
  {
    sleqp_raise(SLEQP_ILLEGAL_ARGUMENT, "Invalid enum option (%d)", option);
  }

  return settings->enum_values[option];
}

SLEQP_RETCODE
sleqp_settings_set_enum_value(SleqpSettings* settings,
                              SLEQP_SETTINGS_ENUM option,
                              int value)
{
  if (!OPT_VALID(option, SLEQP_NUM_ENUM_SETTINGS))
  {
    sleqp_raise(SLEQP_ILLEGAL_ARGUMENT, "Invalid enum option (%d)", option);
  }

  if (!valid_member(option, value))
  {
    sleqp_raise(SLEQP_ILLEGAL_ARGUMENT,
                "Invalid option value (%d) for option %s",
                value,
                sleqp_settings_enum_desc(option));
  }

  settings->enum_values[option] = value;

  return SLEQP_OKAY;
}

int
sleqp_settings_int_value(const SleqpSettings* settings,
                         SLEQP_SETTINGS_INT option)
{
  if (!OPT_VALID(option, SLEQP_NUM_INT_SETTINGS))
  {
    sleqp_raise(SLEQP_ILLEGAL_ARGUMENT, "Invalid int option (%d)", option);
  }

  return settings->int_values[option];
}

SLEQP_RETCODE
sleqp_settings_set_int_value(SleqpSettings* settings,
                             SLEQP_SETTINGS_INT option,
                             int value)
{
  if (!OPT_VALID(option, SLEQP_NUM_INT_SETTINGS))
  {
    sleqp_raise(SLEQP_ILLEGAL_ARGUMENT, "Invalid int option (%d)", option);
  }

  settings->int_values[option] = value;

  return SLEQP_OKAY;
}

bool
sleqp_settings_bool_value(const SleqpSettings* settings,
                          SLEQP_SETTINGS_BOOL option)
{
  if (!OPT_VALID(option, SLEQP_NUM_BOOL_SETTINGS))
  {
    sleqp_raise(SLEQP_ILLEGAL_ARGUMENT, "Invalid bool option (%d)", option);
  }

  return settings->bool_values[option];
}

SLEQP_RETCODE
sleqp_settings_set_bool_value(SleqpSettings* settings,
                              SLEQP_SETTINGS_BOOL option,
                              bool value)
{
  if (!OPT_VALID(option, SLEQP_NUM_BOOL_SETTINGS))
  {
    sleqp_raise(SLEQP_ILLEGAL_ARGUMENT, "Invalid bool option (%d)", option);
  }

  settings->bool_values[option] = value;

  return SLEQP_OKAY;
}

#define BUF_LEN 1024

static void
trim(char* s)
{
  int i;

  while (isspace(*s))
  {
    s++;
  }

  const int size = strlen(s);

  if (size == 0)
  {
    return;
  }

  for (i = size - 1; i >= 0 && (isspace(s[i])); i--)
  {
  }

  s[i + 1] = '\0';
}

void
extract_key_value(char buffer[BUF_LEN],
                  char key[BUF_LEN],
                  char value[BUF_LEN],
                  bool* error,
                  bool* content)
{
  trim(buffer);

  *content = false;
  *error   = false;

  // skip empty line
  if (strlen(buffer) == 0)
  {
    return;
  }

  if (buffer[0] == '#')
  {
    return;
  }

  char* begin = buffer;
  char* end   = buffer;

  for (; *begin && !isspace(*begin); ++begin)
  {
  }

  if (!begin)
  {
    *error = true;
    return;
  }

  assert(isspace(*begin));
  memcpy(key, buffer, (begin - buffer));
  key[begin - buffer] = '\0';

  while (isspace(*begin))
  {
    begin++;
  }

  if (!begin)
  {
    *error = true;
    return;
  }

  end = begin;

  for (; *end && !isspace(*end); ++end)
  {
  }

  memcpy(value, begin, end - begin);
  value[end - begin] = '\0';
  *content           = true;
}

static SLEQP_RETCODE
settings_set_key_value(SleqpSettings* settings,
                       const char* key,
                       const char* value,
                       bool* error)
{
  *error = false;

  // int
  for (int i = 0; i < SLEQP_NUM_INT_SETTINGS; ++i)
  {
    if (!strcmp(key, sleqp_settings_int_name(i)))
    {
      char* end;
      const int v = strtol(value, &end, 10);

      if (*end)
      {
        *error = true;
        return SLEQP_OKAY;
      }

      SLEQP_CALL(sleqp_settings_set_int_value(settings, i, v));
      return SLEQP_OKAY;
    }
  }

  // real
  for (int i = 0; i < SLEQP_NUM_REAL_SETTINGS; ++i)
  {
    if (!strcmp(key, sleqp_settings_real_name(i)))
    {
      char* end;
      const double v = strtod(value, &end);

      if (*end)
      {
        *error = true;
        return SLEQP_OKAY;
      }

      SLEQP_CALL(sleqp_settings_set_real_value(settings, i, v));
      return SLEQP_OKAY;
    }
  }

  // bool
  for (int i = 0; i < SLEQP_NUM_BOOL_SETTINGS; ++i)
  {
    if (!strcmp(key, sleqp_settings_bool_name(i)))
    {
      if (!strcmp(value, "1") || !strcmp(value, "true"))
      {
        SLEQP_CALL(sleqp_settings_set_bool_value(settings, i, true));
        return SLEQP_OKAY;
      }
      else if (!strcmp(value, "0") || !strcmp(value, "false"))
      {
        SLEQP_CALL(sleqp_settings_set_bool_value(settings, i, false));
        return SLEQP_OKAY;
      }
    }
  }

  // enum
  for (int i = 0; i < SLEQP_NUM_ENUM_SETTINGS; ++i)
  {
    if (!strcmp(key, sleqp_settings_enum_name(i)))
    {
      const SleqpEnum* settings_enum = get_enum(i);

      if (!settings_enum)
      {
        *error = true;
        return SLEQP_OKAY;
      }

      const SleqpEnumEntry* entry = settings_enum->entries;

      for (; entry->name; ++entry)
      {
        if (!strcmp(entry->name, value))
        {
          SLEQP_CALL(sleqp_settings_set_enum_value(settings, i, entry->value));
          return SLEQP_OKAY;
        }
      }

      *error = true;
      break;
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_settings_read_file(SleqpSettings* settings, const char* settings_filename)
{
  FILE* file = fopen(settings_filename, "r");

  if (!file)
  {
    sleqp_raise(SLEQP_ILLEGAL_ARGUMENT,
                "Could not open settings file '%s'",
                settings_filename);
  }

  char buffer[BUF_LEN];
  char key[BUF_LEN];
  char value[BUF_LEN];

  int line = 0;

  while (fgets(buffer, BUF_LEN, file))
  {
    ++line;

    bool error;
    bool content;

    extract_key_value(buffer, key, value, &error, &content);

    if (error)
    {
      sleqp_raise(SLEQP_ILLEGAL_ARGUMENT,
                  "Invalid syntax at line %d of settings file '%s'",
                  line,
                  settings_filename);
    }

    if (!content)
    {
      continue;
    }

    SLEQP_CALL(settings_set_key_value(settings, key, value, &error));

    if (error)
    {
      sleqp_raise(
        SLEQP_ILLEGAL_ARGUMENT,
        "Could not set parameter '%s' to '%s' at line %d of settings file '%s'",
        key,
        value,
        line,
        settings_filename);
    }
    else
    {
      sleqp_log_debug("Setting parameter %s to %s", key, value);
    }
  }

  fclose(file);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
settings_free(SleqpSettings** star)
{
  SleqpSettings* settings = *star;

  if (!settings)
  {
    return SLEQP_OKAY;
  }

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_settings_capture(SleqpSettings* settings)
{
  ++settings->refcount;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_settings_release(SleqpSettings** star)
{
  SleqpSettings* settings = *star;

  if (!settings)
  {
    return SLEQP_OKAY;
  }

  if (--settings->refcount == 0)
  {
    SLEQP_CALL(settings_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
