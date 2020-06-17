#include "sleqp_options.h"

#include "sleqp_mem.h"

#define DERIV_CHECK_DEFAULT SLEQP_DERIV_CHECK_SKIP
#define HESSIAN_EVAL_DEFAULT SLEQP_HESSIAN_EVAL_EXACT
#define QUASI_NEWTON_SIZE_DEFAULT 5

struct SleqpOptions
{
  SLEQP_DERIV_CHECK deriv_check;
  SLEQP_HESSIAN_EVAL hessian_eval;
  int quasi_newton_size;
};

SLEQP_RETCODE sleqp_options_create(SleqpOptions** star)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpOptions* options = *star;

  options->deriv_check = DERIV_CHECK_DEFAULT;
  options->hessian_eval = HESSIAN_EVAL_DEFAULT;
  options->quasi_newton_size = QUASI_NEWTON_SIZE_DEFAULT;

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

int sleqp_options_get_quasi_newton_num_iterates(const SleqpOptions* options)
{
  return options->quasi_newton_size;
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
