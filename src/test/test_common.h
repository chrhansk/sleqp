#ifndef TEST_COMMON_H
#define TEST_COMMON_H

#include <check.h>

#include "sleqp_types.h"

#define ASSERT_CALL(x)                          \
  do                                            \
  {                                             \
    SLEQP_RETCODE _retcode_ = (x);              \
    ck_assert_int_eq(_retcode_, SLEQP_OKAY);    \
  }                                             \
  while(0)

#endif /* TEST_COMMON_H */
