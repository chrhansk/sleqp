#ifndef SLEQP_CUTEST_DEFS_H
#define SLEQP_CUTEST_DEFS_H

#include <cutest.h>

#ifdef __cplusplus
extern "C" {
#endif

  extern const logical cutest_true;
  extern const logical cutest_false;

#define SLEQP_CUTEST_CHECK_STATUS(status)       \
  do                                            \
  {                                             \
    if(status)                                  \
    {                                           \
      sleqp_log_error("Error in CUTest call");  \
      return SLEQP_INTERNAL_ERROR;              \
    }                                           \
  }                                             \
  while(0)

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_CUTEST_DEFS_H */
