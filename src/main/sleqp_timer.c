#include "sleqp_timer.h"

#include <time.h>

#include "sleqp_mem.h"

struct SleqpTimer
{
  clock_t start;
};

SLEQP_RETCODE sleqp_timer_create(SleqpTimer** star)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpTimer* timer = *star;

  SLEQP_CALL(sleqp_timer_reset(timer));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_timer_reset(SleqpTimer* timer)
{
  timer->start = clock();

  return SLEQP_OKAY;
}

double sleqp_timer_elapsed(SleqpTimer* timer)
{
  clock_t end = clock();

  return (end - timer->start) / CLOCKS_PER_SEC;
}

SLEQP_RETCODE sleqp_timer_free(SleqpTimer** star)
{
  sleqp_free(star);

  return SLEQP_OKAY;
}
