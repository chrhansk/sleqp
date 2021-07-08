#ifndef SLEQP_TIMER_H
#define SLEQP_TIMER_H

/**
 * @file timer.h
 * @brief Definition of a timer.
 **/

#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpTimer SleqpTimer;

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_timer_create(SleqpTimer** star);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_timer_start(SleqpTimer* timer);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_timer_reset(SleqpTimer* timer);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_timer_stop(SleqpTimer* timer);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_timer_add(SleqpTimer* timer, double value);

  double sleqp_timer_elapsed(SleqpTimer* timer);

  double sleqp_timer_get_avg(SleqpTimer* timer);

  double sleqp_timer_get_ttl(SleqpTimer* timer);

  double sleqp_timer_get_std(SleqpTimer* timer);

  int sleqp_timer_get_num_runs(SleqpTimer* timer);

  SLEQP_NODISCARD SLEQP_RETCODE sleqp_timer_free(SleqpTimer** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_TIMER_H */
