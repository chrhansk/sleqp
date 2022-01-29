#ifndef SLEQP_CUTEST_TYPES_H
#define SLEQP_CUTEST_TYPES_H

#include <cutest.h>

#include "error.h"

extern const logical cutest_true;
extern const logical cutest_false;

extern integer cutest_io_buffer;
extern integer cutest_iout;

#define SLEQP_CUTEST_CHECK_STATUS(status)                                      \
  do                                                                           \
  {                                                                            \
    if (status)                                                                \
    {                                                                          \
      sleqp_raise(SLEQP_FUNC_EVAL_ERROR, "Error in CUTest call");              \
    }                                                                          \
  } while (0)

#endif /* SLEQP_CUTEST_TYPES_H */
