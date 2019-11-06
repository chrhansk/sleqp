#include "sleqp_options.h"

#include "sleqp_mem.h"

struct SleqpOptions
{
  SLEQP_DERIV_CHECK deriv_check;
  SLEQP_HESSIAN_EVAL hessian_eval;
};

SLEQP_RETCODE sleqp_options_create(SleqpOptions** star)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpOptions* options = *star;

  options->deriv_check = SLEQP_DERIV_CHECK_SKIP;
  options->hessian_eval = SLEQP_HESSIAN_EVAL_EXACT;

  return SLEQP_OKAY;
}


SLEQP_DERIV_CHECK sleqp_options_get_deriv_check(const SleqpOptions* options)
{
  return options->deriv_check;
}

SLEQP_HESSIAN_EVAL sleqp_options_get_hessian_eval(const SleqpOptions* options)
{
  return options->hessian_eval;
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
