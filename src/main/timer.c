#include "timer.h"

#include <math.h>
#include <time.h>

#include "cmp.h"
#include "fail.h"
#include "log.h"
#include "mem.h"

#define BUF_SIZE 512

struct SleqpTimer
{
  clock_t start;
  int num_runs;
  int running;

  double total_elapsed;
  double total_elapsed_squared;

  double last_elapsed;
};

typedef enum
{
  Day = 0,
  Hour,
  Minute,
  Second,
  Milli,
  Micro,
  NumUnits
} TimeUnit;

typedef struct
{
  const double factor;
  const char* suffix;
  int len;
} TimeConversion;

const int suffix_format_len = 3;

static const TimeConversion conversions[NumUnits]
  = {[Day]    = {.factor = 1 / 24., .suffix = "d", .len = 1},
     [Hour]   = {.factor = 1 / 60., .suffix = "h", .len = 1},
     [Minute] = {.factor = 1 / 60., .suffix = "min", .len = 3},
     [Second] = {.factor = 1., .suffix = "s", .len = 1},
     [Milli]  = {.factor = 1000., .suffix = "ms", .len = 2},
     [Micro]  = {.factor = 1000., .suffix = SLEQP_SYMBOL_MU "s", .len = 2}};

SLEQP_RETCODE
sleqp_timer_create(SleqpTimer** star)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpTimer* timer = *star;

  SLEQP_CALL(sleqp_timer_reset(timer));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_timer_start(SleqpTimer* timer)
{
  assert(!timer->running);

  timer->running = true;
  timer->start   = clock();

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_timer_reset(SleqpTimer* timer)
{
  *timer = (SleqpTimer){0};

  return SLEQP_OKAY;
}

static double
current_elapsed(SleqpTimer* timer)
{
  if (!timer->running)
  {
    return 0.;
  }

  clock_t end = clock();

  double elapsed = ((double)(end - timer->start)) / CLOCKS_PER_SEC;

  assert(elapsed >= 0.);

  return elapsed;
}

SLEQP_RETCODE
sleqp_timer_stop(SleqpTimer* timer)
{
  assert(timer->running);

  const double elapsed = current_elapsed(timer);

  ++timer->num_runs;

  timer->last_elapsed = elapsed;

  timer->total_elapsed += elapsed;
  timer->total_elapsed_squared += elapsed * elapsed;

  timer->running = false;

  return SLEQP_OKAY;
}

static double
currently_elapsed(SleqpTimer* timer)
{
  if (!timer->running)
  {
    return 0.;
  }

  clock_t end = clock();

  return ((double)(end - timer->start)) / CLOCKS_PER_SEC;
}

SLEQP_RETCODE
sleqp_timer_add(SleqpTimer* timer, double value)
{
  timer->last_elapsed += value;

  timer->total_elapsed += value;
  timer->total_elapsed_squared += value * value;

  return SLEQP_OKAY;
}

double
sleqp_timer_elapsed(SleqpTimer* timer)
{
  return timer->last_elapsed + currently_elapsed(timer);
}

double
sleqp_timer_get_avg(SleqpTimer* timer)
{
  assert(!timer->running);

  if (timer->num_runs == 0)
  {
    return 0.;
  }

  return timer->total_elapsed / timer->num_runs;
}

double
sleqp_timer_get_ttl(SleqpTimer* timer)
{
  return timer->total_elapsed + currently_elapsed(timer);
}

double
sleqp_timer_get_std(SleqpTimer* timer)
{
  if (timer->num_runs == 0)
  {
    return 0.;
  }

  double avg = sleqp_timer_get_avg(timer);

  double var = timer->total_elapsed_squared / timer->num_runs - avg * avg;

  return sqrt(var);
}

int
sleqp_timer_get_num_runs(SleqpTimer* timer)
{
  return timer->num_runs;
}

static SLEQP_RETCODE
round_up_time(double* time, const char** suffix, int* len)
{
  for (int idx = Second; idx >= 0; --idx)
  {
    const TimeConversion* conversion = conversions + idx;

    const double next_time = (*time) * conversion->factor;

    if (next_time >= 1.)
    {
      *time   = next_time;
      *suffix = conversion->suffix;
      *len    = conversion->len;
    }
    else
    {
      break;
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
round_down_time(double* time, const char** suffix, int* len)
{
  for (int idx = Second; idx < NumUnits; ++idx)
  {
    if (*time >= 1.)
    {
      break;
    }

    const TimeConversion* conversion = conversions + idx;

    *time *= conversion->factor;
    *suffix = conversion->suffix;
    *len    = conversion->len;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
round_time(double* time, const char** suffix, int* len)
{
  if (*time == 0.)
  {
    *suffix = conversions[Second].suffix;
    return SLEQP_OKAY;
  }

  if (*time >= 1.)
  {
    return round_up_time(time, suffix, len);
  }
  else
  {
    return round_down_time(time, suffix, len);
  }
}

SLEQP_RETCODE
sleqp_timer_display(SleqpTimer* timer,
                    const char* description,
                    double total_elapsed)
{
  char buffer[BUF_SIZE];

  const int num_runs = sleqp_timer_get_num_runs(timer);

  double avg_time   = sleqp_timer_get_avg(timer);
  double total_time = sleqp_timer_get_ttl(timer);

  if (num_runs == 0)
  {
    snprintf(buffer, BUF_SIZE, "%40s: %5d", description, num_runs);
  }
  else
  {
    const double percent = (total_time / total_elapsed) * 100.;

    const char* total_suffix = NULL;
    int total_len;
    SLEQP_CALL(round_time(&total_time, &total_suffix, &total_len));

    if (num_runs == 1)
    {
      snprintf(buffer,
               BUF_SIZE,
               "%40s: %5d (%6.2f %*s%s = %6.2f%%)",
               description,
               num_runs,
               total_time,
               suffix_format_len - total_len,
               "",
               total_suffix,
               percent);
    }
    else
    {
      const char* avg_suffix = NULL;
      int avg_len;
      SLEQP_CALL(round_time(&avg_time, &avg_suffix, &avg_len));

      snprintf(buffer,
               BUF_SIZE,
               "%40s: %5d (%6.2f %*s%s avg, %8.2f %*s%s total = %6.2f%%)",
               description,
               num_runs,
               avg_time,
               suffix_format_len - avg_len,
               "",
               avg_suffix,
               total_time,
               suffix_format_len - total_len,
               "",
               total_suffix,
               percent);
    }
  }

  sleqp_log_info("%s", buffer);

  return SLEQP_OKAY;
}

double
sleqp_timer_remaining_time(SleqpTimer* timer, double time_limit)
{
  return sleqp_remaining_time(sleqp_timer_get_ttl(timer), time_limit);
}

double
sleqp_remaining_time(double elapsed_time, double time_limit)
{
  if (time_limit != SLEQP_NONE)
  {
    double remaining_time = time_limit - elapsed_time;

    return SLEQP_MAX(remaining_time, 0.);
  }

  return SLEQP_NONE;
}

bool
sleqp_timer_exhausted_time_limit(SleqpTimer* timer, double time_limit)
{
  if (time_limit == SLEQP_NONE)
  {
    return false;
  }

  return sleqp_timer_remaining_time(timer, time_limit) <= 0.;
}

SLEQP_RETCODE
sleqp_timer_free(SleqpTimer** star)
{
  sleqp_free(star);

  return SLEQP_OKAY;
}
