#ifndef SLEQP_CUTEST_DRIVER_H
#define SLEQP_CUTEST_DRIVER_H

#include "sleqp_cutest_options.h"

#ifdef __cplusplus
extern "C" {
#endif

  int sleqp_cutest_run(const char* filename,
                       const char* probname,
                       const SleqpCutestOptions* cutest_options);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_CUTEST_DRIVER_H */
