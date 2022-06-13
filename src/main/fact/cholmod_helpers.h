#ifndef SLEQP_CHOLMOD_HELPERS_H
#define SLEQP_CHOLMOD_HELPERS_H

#include "types.h"

void
sleqp_cholmod_report_error(int status,
                           const char* file,
                           int line,
                           const char* message);

SLEQP_RETCODE
sleqp_cholmod_error_string(int status, const char** message);

#define SLEQP_CHOLMOD_ERROR_CHECK(common)                                      \
  do                                                                           \
  {                                                                            \
    if ((common)->status < 0)                                                  \
    {                                                                          \
      const char* message;                                                     \
      SLEQP_CALL(sleqp_cholmod_error_string((common)->status, &message));      \
      sleqp_raise(SLEQP_INTERNAL_ERROR,                                        \
                  "Caught CHOLMOD error <%d> (%s)",                            \
                  (common)->status,                                            \
                  message);                                                    \
    }                                                                          \
  } while (false)

#endif /* SLEQP_CHOLMOD_HELPERS_H */
