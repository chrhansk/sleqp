#include "sleqp_cutest_options.h"

#define TIME_LIMIT_DEFAULT 3600.
#define MAX_NUM_THREADS_DEFAULT SLEQP_NONE

SLEQP_RETCODE sleqp_cutest_options_default(SleqpCutestOptions* cutest_options)
{
  *cutest_options = (SleqpCutestOptions) {0};

  cutest_options->time_limit = TIME_LIMIT_DEFAULT;

  cutest_options->max_num_threads = MAX_NUM_THREADS_DEFAULT;

  return SLEQP_OKAY;
}
