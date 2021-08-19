#ifndef SLEQP_CUTEST_OPTIONS_H
#define SLEQP_CUTEST_OPTIONS_H

#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  bool enable_logging;
  bool enable_preprocessing;
  bool force_nonlinear_constraints;
  int max_num_threads;
  const char* output;

} SleqpCutestOptions;

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_CUTEST_OPTIONS_H */
