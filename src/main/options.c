#include "options.h"

#include <assert.h>
#include <fenv.h>
#include <math.h>

#include "enum.h"
#include "log.h"
#include "mem.h"
#include "types.h"

#define PERFORM_NEWTON_DEFAULT true
#define PERFORM_SOC_DEFAULT true
#define USE_QUADRATIC_MODEL_DEFAULT true
#define ALWAYS_WARM_START_LP_DEFAULT true
#define ENABLE_PREPROCESSOR_DEFAULT false
#define ENABLE_RESTORATION_PHASE_DEFAULT false

#define PARAMETRIC_CAUCHY_DEFAULT SLEQP_PARAMETRIC_CAUCHY_DISABLED
#define INITIAL_TR_CHOICE_DEFAULT SLEQP_INITIAL_TR_CHOICE_NARROW
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

struct SleqpOptions
{
  int refcount;

  int enum_values[SLEQP_NUM_ENUM_OPTIONS];
  int int_values[SLEQP_NUM_INT_OPTIONS];
  bool bool_values[SLEQP_NUM_BOOL_OPTIONS];

  SLEQP_BFGS_SIZING bfgs_sizing;
  SLEQP_TR_SOLVER tr_solver;
};

SLEQP_RETCODE
sleqp_options_create(SleqpOptions** star)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpOptions* options = *star;

  *options = (SleqpOptions){0};

  *options = (SleqpOptions){
    .refcount = 1,
    .enum_values
    = {[SLEQP_OPTION_ENUM_DERIV_CHECK]          = DERIV_CHECK_DEFAULT,
       [SLEQP_OPTION_ENUM_HESS_EVAL]            = HESS_EVAL_DEFAULT,
       [SLEQP_OPTION_ENUM_DUAL_ESTIMATION_TYPE] = DUAL_ESTIMATION_TYPE_DEFAULT,
       [SLEQP_OPTION_ENUM_FLOAT_WARNING_FLAGS]  = FLOAT_WARN_FLAGS_DEFAULT,
       [SLEQP_OPTION_ENUM_FLOAT_ERROR_FLAGS]    = FLOAT_ERR_FLAGS_DEFAULT,
       [SLEQP_OPTION_ENUM_BFGS_SIZING]          = BFGS_SIZING_DEFAULT,
       [SLEQP_OPTION_ENUM_TR_SOLVER]            = TR_SOLVER_DEFAULT,
       [SLEQP_OPTION_ENUM_POLISHING_TYPE]       = POLISHING_TYPE_DEFAULT,
       [SLEQP_OPTION_ENUM_STEP_RULE]            = STEP_RULE_DEFAULT,
       [SLEQP_OPTION_ENUM_LINESEARCH]           = LINESEARCH_DEFAULT,
       [SLEQP_OPTION_ENUM_PARAMETRIC_CAUCHY]    = PARAMETRIC_CAUCHY_DEFAULT,
       [SLEQP_OPTION_ENUM_INITIAL_TR_CHOICE]    = INITIAL_TR_CHOICE_DEFAULT},
    .int_values
    = {[SLEQP_OPTION_INT_NUM_QUASI_NEWTON_ITERATES] = QUASI_NEWTON_SIZE_DEFAULT,
       [SLEQP_OPTION_INT_MAX_NEWTON_ITERATIONS] = MAX_NEWTON_ITERATIONS_DEFAULT,
       [SLEQP_OPTION_INT_NUM_THREADS]           = NUM_THREADS_DEFAULT},
    .bool_values
    = {[SLEQP_OPTION_BOOL_PERFORM_NEWTON_STEP]  = PERFORM_NEWTON_DEFAULT,
       [SLEQP_OPTION_BOOL_PERFORM_SOC]          = PERFORM_SOC_DEFAULT,
       [SLEQP_OPTION_BOOL_USE_QUADRATIC_MODEL]  = USE_QUADRATIC_MODEL_DEFAULT,
       [SLEQP_OPTION_BOOL_ALWAYS_WARM_START_LP] = ALWAYS_WARM_START_LP_DEFAULT,
       [SLEQP_OPTION_BOOL_ENABLE_RESTORATION_PHASE]
       = ENABLE_RESTORATION_PHASE_DEFAULT,
       [SLEQP_OPTION_BOOL_ENABLE_PREPROCESSOR] = ENABLE_PREPROCESSOR_DEFAULT}};

  return SLEQP_OKAY;
}

static bool
valid_member(SLEQP_OPTION_ENUM option, int value)
{
  switch (option)
  {
  case SLEQP_OPTION_ENUM_DERIV_CHECK:
    return sleqp_enum_member(sleqp_enum_deriv_check(), value);
  case SLEQP_OPTION_ENUM_HESS_EVAL:
    return sleqp_enum_member(sleqp_enum_hess_eval(), value);
  case SLEQP_OPTION_ENUM_DUAL_ESTIMATION_TYPE:
    return sleqp_enum_member(sleqp_enum_dual_estimation(), value);
  case SLEQP_OPTION_ENUM_FLOAT_WARNING_FLAGS:
  case SLEQP_OPTION_ENUM_FLOAT_ERROR_FLAGS:
    return true;
  case SLEQP_OPTION_ENUM_BFGS_SIZING:
    return sleqp_enum_member(sleqp_enum_bfgs_sizing(), value);
  case SLEQP_OPTION_ENUM_TR_SOLVER:
    return sleqp_enum_member(sleqp_enum_tr_solver(), value);
  case SLEQP_OPTION_ENUM_POLISHING_TYPE:
    return sleqp_enum_member(sleqp_enum_polishing_type(), value);
  case SLEQP_OPTION_ENUM_STEP_RULE:
    return sleqp_enum_member(sleqp_enum_step_rule(), value);
  case SLEQP_OPTION_ENUM_LINESEARCH:
    return sleqp_enum_member(sleqp_enum_linesearch(), value);
  case SLEQP_OPTION_ENUM_PARAMETRIC_CAUCHY:
    return sleqp_enum_member(sleqp_enum_parametric_cauchy(), value);
  case SLEQP_OPTION_ENUM_INITIAL_TR_CHOICE:
    return sleqp_enum_member(sleqp_enum_initial_tr(), value);
  default:
    assert(0);
  }
  return false;
}

int
sleqp_options_enum_value(const SleqpOptions* options, SLEQP_OPTION_ENUM option)
{
  assert(option >= 0);
  assert(option < SLEQP_NUM_ENUM_OPTIONS);

  return options->enum_values[option];
}

SLEQP_RETCODE
sleqp_options_set_enum_value(SleqpOptions* options,
                             SLEQP_OPTION_ENUM option,
                             int value)
{
  assert(option >= 0);
  assert(option < SLEQP_NUM_ENUM_OPTIONS);

  if (!valid_member(option, value))
  {
    return SLEQP_ILLEGAL_ARGUMENT;
  }

  options->enum_values[option] = value;

  return SLEQP_OKAY;
}

int
sleqp_options_int_value(const SleqpOptions* options, SLEQP_OPTION_INT option)
{
  assert(option >= 0);
  assert(option < SLEQP_NUM_INT_OPTIONS);

  return options->int_values[option];
}

SLEQP_RETCODE
sleqp_options_set_int_value(SleqpOptions* options,
                            SLEQP_OPTION_INT option,
                            int value)
{
  assert(option >= 0);
  assert(option < SLEQP_NUM_INT_OPTIONS);

  options->int_values[option] = value;

  return SLEQP_OKAY;
}

bool
sleqp_options_bool_value(const SleqpOptions* options, SLEQP_OPTION_BOOL option)
{
  assert(option >= 0);
  assert(option < SLEQP_NUM_BOOL_OPTIONS);

  return options->bool_values[option];
}

SLEQP_RETCODE
sleqp_options_set_bool_value(SleqpOptions* options,
                             SLEQP_OPTION_BOOL option,
                             bool value)
{
  assert(option >= 0);
  assert(option < SLEQP_NUM_BOOL_OPTIONS);

  options->bool_values[option] = value;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
options_free(SleqpOptions** star)
{
  SleqpOptions* options = *star;

  if (!options)
  {
    return SLEQP_OKAY;
  }

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_options_capture(SleqpOptions* options)
{
  ++options->refcount;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_options_release(SleqpOptions** star)
{
  SleqpOptions* options = *star;

  if (!options)
  {
    return SLEQP_OKAY;
  }

  if (--options->refcount == 0)
  {
    SLEQP_CALL(options_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
