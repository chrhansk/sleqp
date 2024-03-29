#ifndef SLEQP_TIMER_H
#define SLEQP_TIMER_H

/**
 * @file timer.h
 * @brief Definition of a timer.
 **/

#include "types.h"

typedef struct SleqpTimer SleqpTimer;

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_timer_create(SleqpTimer** star);

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_timer_start(SleqpTimer* timer);

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_timer_reset(SleqpTimer* timer);

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_timer_stop(SleqpTimer* timer);

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_timer_add(SleqpTimer* timer, double value);

double
sleqp_timer_elapsed(SleqpTimer* timer);

double
sleqp_timer_get_avg(SleqpTimer* timer);

double
sleqp_timer_get_ttl(SleqpTimer* timer);

double
sleqp_timer_get_std(SleqpTimer* timer);

int
sleqp_timer_get_num_runs(SleqpTimer* timer);

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_timer_display(SleqpTimer* timer,
                    const char* description,
                    double total_elapsed);

double
sleqp_timer_remaining_time(SleqpTimer* timer, double time_limit);

double
sleqp_remaining_time(double elapsed_time, double time_limit);

bool
sleqp_timer_exhausted_time_limit(SleqpTimer* timer, double time_limit);

SLEQP_WARNUNUSED SLEQP_RETCODE
sleqp_timer_free(SleqpTimer** star);

#endif /* SLEQP_TIMER_H */
