#include "sleqp_options.h"

#include <fenv.h>
#include <math.h>

#include "sleqp_log.h"
#include "sleqp_mem.h"

#define PERFORM_NEWTON_DEFAULT             true
#define PERFORM_SOC_DEFAULT                true
#define USE_QUADRATIC_MODEL_DEFAULT        true
#define DERIV_CHECK_DEFAULT                SLEQP_DERIV_CHECK_SKIP
#define HESSIAN_EVAL_DEFAULT               SLEQP_HESSIAN_EVAL_EXACT
#define SLEQP_DUAL_ESTIMATION_TYPE_DEFAULT SLEQP_DUAL_ESTIMATION_TYPE_LSQ
#define QUASI_NEWTON_SIZE_DEFAULT          5
#define MAX_NEWTON_ITERATIONS_DEFAULT      100
#define FLOAT_WARN_FLAGS_DEFAULT           FE_ALL_EXCEPT
#define FLOAT_ERR_FLAGS_DEFAULT            (FE_OVERFLOW | FE_DIVBYZERO | FE_INVALID)
#define BFGS_SIZING_DEFAULT                SLEQP_BFGS_SIZING_CENTERED_OL

#define CHECK_FLOAT_ENV                                                                             \
  do                                                                                                \
  {                                                                                                 \
    if(!(math_errhandling & MATH_ERREXCEPT))                                                        \
    {                                                                                               \
      sleqp_log_warn("Float point error handling is not supported, setting options has no effect"); \
    }                                                                                               \
  }                                                                                                 \
  while(false)

struct SleqpOptions
{
  int refcount;

  bool perform_newton_step;
  bool perform_soc;
  bool use_quadratic_model;
  SLEQP_DERIV_CHECK deriv_check;
  SLEQP_HESSIAN_EVAL hessian_eval;
  SLEQP_DUAL_ESTIMATION_TYPE dual_estimation_type;
  int quasi_newton_size;
  int max_newton_iterations;

  int float_warn_flags;
  int float_err_flags;

  SLEQP_BFGS_SIZING bfgs_sizing;
};

SLEQP_RETCODE sleqp_options_create(SleqpOptions** star)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpOptions* options = *star;

  *options = (SleqpOptions){0};

  *options = (SleqpOptions) {
    .refcount = 1,

    .perform_newton_step = PERFORM_NEWTON_DEFAULT,
    .perform_soc = PERFORM_SOC_DEFAULT,
    .use_quadratic_model = USE_QUADRATIC_MODEL_DEFAULT,
    .deriv_check = DERIV_CHECK_DEFAULT,
    .hessian_eval = HESSIAN_EVAL_DEFAULT,
    .dual_estimation_type = SLEQP_DUAL_ESTIMATION_TYPE_DEFAULT,
    .quasi_newton_size = QUASI_NEWTON_SIZE_DEFAULT,
    .max_newton_iterations = MAX_NEWTON_ITERATIONS_DEFAULT,

    .float_warn_flags = FLOAT_WARN_FLAGS_DEFAULT,
    .float_err_flags = FLOAT_ERR_FLAGS_DEFAULT,

    .bfgs_sizing = BFGS_SIZING_DEFAULT,
  };

  return SLEQP_OKAY;
}

bool sleqp_options_get_perform_newton_step(const SleqpOptions* options)
{
  return options->perform_newton_step;
}

bool sleqp_options_get_perform_soc(const SleqpOptions* options)
{
  return options->perform_soc;
}

bool sleqp_options_get_use_quadratic_model(const SleqpOptions* options)
{
  return options->use_quadratic_model;
}


SLEQP_DERIV_CHECK sleqp_options_get_deriv_check(const SleqpOptions* options)
{
  return options->deriv_check;
}

SLEQP_HESSIAN_EVAL sleqp_options_get_hessian_eval(const SleqpOptions* options)
{
  return options->hessian_eval;
}

SLEQP_DUAL_ESTIMATION_TYPE
sleqp_options_get_dual_estimation_type(const SleqpOptions* options)
{
  return options->dual_estimation_type;
}

int sleqp_options_get_quasi_newton_num_iterates(const SleqpOptions* options)
{
  return options->quasi_newton_size;
}

int sleqp_options_get_max_newton_iterations(const SleqpOptions* options)
{
  return options->max_newton_iterations;
}

int sleqp_options_get_float_warning_flags(const SleqpOptions* options)
{
  return options->float_warn_flags;
}

int sleqp_options_get_float_error_flags(const SleqpOptions* options)
{
  return options->float_err_flags;
}

SLEQP_BFGS_SIZING sleqp_options_get_bfgs_sizing(const SleqpOptions* options)
{
  return options->bfgs_sizing;
}

SLEQP_RETCODE sleqp_options_set_perform_newton_step(SleqpOptions* options,
                                                    bool value)
{
  options->perform_newton_step = value;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_options_set_perform_soc(SleqpOptions* options,
                                            bool value)
{
  options->perform_soc = value;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_options_set_use_quadratic_model(SleqpOptions* options, bool value)
{
  options->use_quadratic_model = value;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_options_set_deriv_check(SleqpOptions* options,
                                            SLEQP_DERIV_CHECK value)
{
  options->deriv_check = value;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_options_set_hessian_eval(SleqpOptions* options,
                                             SLEQP_HESSIAN_EVAL value)
{
  options->hessian_eval = value;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_options_set_dual_estimation_type(SleqpOptions* options,
                                                     SLEQP_DUAL_ESTIMATION_TYPE dual_estimation_type)
{
  options->dual_estimation_type = dual_estimation_type;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_options_set_quasi_newton_num_iterates(SleqpOptions* options,
                                                          int size)
{
  if(size <= 0)
  {
    return SLEQP_ILLEGAL_ARGUMENT;
  }

  options->quasi_newton_size = size;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_options_set_max_newton_iterations(SleqpOptions* options, int iterations)
{
  if((iterations < 0) && (iterations != -1))
  {
    return SLEQP_ILLEGAL_ARGUMENT;
  }

  options->max_newton_iterations = iterations;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_options_set_float_warning_flags(SleqpOptions* options, int flags)
{
  CHECK_FLOAT_ENV;

  options->float_warn_flags = flags;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_options_set_float_error_flags(SleqpOptions* options, int flags)
{
  CHECK_FLOAT_ENV;

  options->float_err_flags = flags;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_options_set_bfgs_sizing(SleqpOptions* options, SLEQP_BFGS_SIZING sizing)
{
  options->bfgs_sizing = sizing;

  return SLEQP_OKAY;
}

SLEQP_RETCODE options_free(SleqpOptions** star)
{
  SleqpOptions* options = *star;

  if(!options)
  {
    return SLEQP_OKAY;
  }

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_options_capture(SleqpOptions* options)
{
  ++options->refcount;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_options_release(SleqpOptions** star)
{
  SleqpOptions* options = *star;

  if(!options)
  {
    return SLEQP_OKAY;
  }

  if(--options->refcount == 0)
  {
    SLEQP_CALL(options_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
