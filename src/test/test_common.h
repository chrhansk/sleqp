#ifndef TEST_COMMON_H
#define TEST_COMMON_H

#include "sleqp_types.h"

#define CK_ASSERT_SLEQP_CALL(x)                 \
  do                                            \
  {                                             \
    SLEQP_RETCODE _retcode_ = (x);              \
    ck_assert_int_eq(_retcode_, SLEQP_OKAY);    \
  }                                             \
  while(0)

#endif /* TEST_COMMON_H */
