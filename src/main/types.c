#include "types.h"

const char*
sleqp_retcode_str(SLEQP_RETCODE retcode)
{
  switch (retcode)
  {
  case SLEQP_NOMEM:
    return "SLEQP_NOMEM";
  case SLEQP_FAILED_ASSERTION:
    return "SLEQP_FAILED_ASSERTION";
  case SLEQP_ILLEGAL_ARGUMENT:
    return "SLEQP_ILLEGAL_ARGUMENT";
  case SLEQP_INVALID_DERIV:
    return "SLEQP_INVALID_DERIV";
  case SLEQP_INTERNAL_ERROR:
    return "SLEQP_INTERNAL_ERROR";
  case SLEQP_MATH_ERROR:
    return "SLEQP_MATH_ERROR";
  case SLEQP_OKAY:
    return "SLEQP_OKAY";
  case SLEQP_ABORT_TIME:
    return "SLEQP_ABORT_TIME";
  }

  return "<unknown>";
}
