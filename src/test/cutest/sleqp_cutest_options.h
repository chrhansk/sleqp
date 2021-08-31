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
    bool no_fork;
    int max_num_threads;
    double time_limit;
    const char* output;

  } SleqpCutestOptions;

  SLEQP_RETCODE sleqp_cutest_options_default(SleqpCutestOptions* cutest_options);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_CUTEST_OPTIONS_H */
