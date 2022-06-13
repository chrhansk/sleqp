#include "cholmod_helpers.h"

#include <cholmod.h>

void
sleqp_cholmod_report_error(int status,
                           const char* file,
                           int line,
                           const char* message)
{
  if (status < 0)
  {
    sleqp_log_error("CHOLMOD: '%s'", message);
  }
  else
  {
    sleqp_log_warn("CHOLMOD: '%s'", message);
  }
}

SLEQP_RETCODE
sleqp_cholmod_error_string(int status, const char** message)
{
  switch (status)
  {
  case CHOLMOD_NOT_INSTALLED:
    (*message) = "method not installed";
    break;
  case CHOLMOD_OUT_OF_MEMORY:
    (*message) = "out of memory";
    break;
  case CHOLMOD_TOO_LARGE:
    (*message) = "integer overflow occured";
    break;
  case CHOLMOD_INVALID:
    (*message) = "invalid input";
    break;
  case CHOLMOD_GPU_PROBLEM:
    (*message) = "GPU fatal error";
    break;
  default:
    (*message) = "unknown error";
  }

  return SLEQP_OKAY;
}
