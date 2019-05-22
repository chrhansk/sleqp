#ifndef SLEQP_TIMER_H
#define SLEQP_TIMER_H

/**
 * @file sleqp_timer.h
 * @brief Definition of a timer.
 **/

#include "sleqp_types.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpTimer SleqpTimer;

  SLEQP_RETCODE sleqp_timer_create(SleqpTimer** star);

  SLEQP_RETCODE sleqp_timer_start(SleqpTimer* timer);

  SLEQP_RETCODE sleqp_timer_reset(SleqpTimer* timer);

  SLEQP_RETCODE sleqp_timer_stop(SleqpTimer* timer);

  double sleqp_timer_elapsed(SleqpTimer* timer);

  double sleqp_timer_get_avg(SleqpTimer* timer);

  double sleqp_timer_get_ttl(SleqpTimer* timer);

  double sleqp_timer_get_std(SleqpTimer* timer);

  int sleqp_timer_get_num_runs(SleqpTimer* timer);

  SLEQP_RETCODE sleqp_timer_free(SleqpTimer** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_TIMER_H */
