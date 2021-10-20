#include "timer.h"

#include <assert.h>
#include <math.h>
#include <time.h>

#include "log.h"
#include "mem.h"

struct SleqpTimer
{
  clock_t start;
  int num_runs;
  int running;

  double total_elapsed;
  double total_elapsed_squared;

  double last_elapsed;
};

SLEQP_RETCODE sleqp_timer_create(SleqpTimer** star)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpTimer* timer = *star;

  SLEQP_CALL(sleqp_timer_reset(timer));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_timer_start(SleqpTimer* timer)
{
  assert(!timer->running);

  timer->running = true;
  timer->start = clock();

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_timer_reset(SleqpTimer* timer)
{
  *timer = (SleqpTimer){0};

  return SLEQP_OKAY;
}

static
double current_elapsed(SleqpTimer* timer)
{
  if(!timer->running)
  {
    return 0.;
  }

  clock_t end = clock();

  double elapsed = ((double) (end - timer->start)) / CLOCKS_PER_SEC;

  assert(elapsed >= 0.);

  return elapsed;
}

SLEQP_RETCODE sleqp_timer_stop(SleqpTimer* timer)
{
  assert(timer->running);

  const double elapsed = current_elapsed(timer);

  ++timer->num_runs;

  timer->last_elapsed = elapsed;

  timer->total_elapsed += elapsed;
  timer->total_elapsed_squared += elapsed*elapsed;

  timer->running = false;

  return SLEQP_OKAY;
}

static double
currently_elapsed(SleqpTimer* timer)
{
  if(!timer->running)
  {
    return 0.;
  }

  clock_t end = clock();

  return ((double) (end - timer->start)) / CLOCKS_PER_SEC;
}

SLEQP_RETCODE sleqp_timer_add(SleqpTimer* timer, double value)
{
  timer->last_elapsed += value;

  timer->total_elapsed += value;
  timer->total_elapsed_squared += value*value;

  return SLEQP_OKAY;
}

double sleqp_timer_elapsed(SleqpTimer* timer)
{
  return timer->last_elapsed + currently_elapsed(timer);
}

double sleqp_timer_get_avg(SleqpTimer* timer)
{
  assert(!timer->running);

  if(timer->num_runs == 0)
  {
    return 0.;
  }

  return timer->total_elapsed / timer->num_runs;
}

double sleqp_timer_get_ttl(SleqpTimer* timer)
{
  return timer->total_elapsed + currently_elapsed(timer);
}

double sleqp_timer_get_std(SleqpTimer* timer)
{
  if(timer->num_runs == 0)
  {
    return 0.;
  }

  double avg = sleqp_timer_get_avg(timer);

  double var = timer->total_elapsed_squared / timer->num_runs - avg*avg;

  return sqrt(var);
}

int sleqp_timer_get_num_runs(SleqpTimer* timer)
{
  return timer->num_runs;
}

SLEQP_RETCODE sleqp_timer_free(SleqpTimer** star)
{
  sleqp_free(star);

  return SLEQP_OKAY;
}
