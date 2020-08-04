#include "sleqp_options.h"

#include "sleqp_mem.h"

#define PERFORM_NEWTON_DEFAULT             true
#define PERFORM_SOC_DEFAULT                true
#define DERIV_CHECK_DEFAULT                SLEQP_DERIV_CHECK_SKIP
#define HESSIAN_EVAL_DEFAULT               SLEQP_HESSIAN_EVAL_EXACT
#define SLEQP_DUAL_ESTIMATION_TYPE_DEFAULT SLEQP_DUAL_ESTIMATION_TYPE_LSQ
#define QUASI_NEWTON_SIZE_DEFAULT          5

struct SleqpOptions
{
  bool perform_newton_step;
  bool perform_soc;
  SLEQP_DERIV_CHECK deriv_check;
  SLEQP_HESSIAN_EVAL hessian_eval;
  SLEQP_DUAL_ESTIMATION_TYPE dual_estimation_type;
  int quasi_newton_size;
};

SLEQP_RETCODE sleqp_options_create(SleqpOptions** star)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpOptions* options = *star;

  options->perform_newton_step = PERFORM_NEWTON_DEFAULT;
  options->perform_soc = PERFORM_SOC_DEFAULT;
  options->deriv_check = DERIV_CHECK_DEFAULT;
  options->hessian_eval = HESSIAN_EVAL_DEFAULT;
  options->dual_estimation_type = SLEQP_DUAL_ESTIMATION_TYPE_DEFAULT;
  options->quasi_newton_size = QUASI_NEWTON_SIZE_DEFAULT;

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

SLEQP_RETCODE sleqp_options_free(SleqpOptions** star)
{
  SleqpOptions* options = *star;

  if(!options)
  {
    return SLEQP_OKAY;
  }

  sleqp_free(star);

  return SLEQP_OKAY;
}
