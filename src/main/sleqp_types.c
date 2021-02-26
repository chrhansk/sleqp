#include "sleqp_types.h"

const char* const sleqp_retcode_names[] = {
  [SLEQP_OKAY]             = "SLEQP_OKAY",
  [SLEQP_NOMEM]            = "SLEQP_NOMEM",
  [SLEQP_FAILED_ASSERTION] = "SLEQP_FAILED_ASSERTION",
  [SLEQP_ILLEGAL_ARGUMENT] = "SLEQP_ILLEGAL_ARGUMENT",
  [SLEQP_INVALID_DERIV]    = "SLEQP_INVALID_DERIV",
  [SLEQP_INTERNAL_ERROR]   = "SLEQP_INTERNAL_ERROR",
  [SLEQP_MATH_ERROR]       = "SLEQP_MATH_ERROR"
};
