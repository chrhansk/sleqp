#include "log.h"

#include <pthread.h>

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TIME_BUF_SIZE 128
#define TOTAL_BUF_SIZE 2048
#define EXTENDED_BUF_SIZE 4096

struct LevelInfo
{
  const char* name;
  const char* color;
};

static struct LevelInfo const level_infos[SLEQP_NUM_LOG_LEVELS]
  = {{"silent", SLEQP_FORMAT_RED},
     {"error", SLEQP_FORMAT_RED},
     {"warn", SLEQP_FORMAT_YELLOW},
     {"info", SLEQP_FORMAT_GREEN},
     {"debug", SLEQP_FORMAT_BLUE}};

#if SLEQP_DEBUG
static SLEQP_LOG_LEVEL level = SLEQP_LOG_INFO;
#else
static SLEQP_LOG_LEVEL level = SLEQP_LOG_DEBUG;
#endif

SLEQP_LOG_LEVEL
sleqp_log_level()
{
  return level;
}

void
sleqp_log_set_level(SLEQP_LOG_LEVEL value)
{
  level = value;
}

struct tm*
localtime_synced(const time_t* timer, struct tm* result)
{
  static pthread_mutex_t localtime_mutex = PTHREAD_MUTEX_INITIALIZER;

  struct tm* r = NULL;

  if (pthread_mutex_lock(&localtime_mutex) == 0)
  {
    r = localtime(timer);
    if (r != 0)
    {
      *result = *r;
      r       = result;
    }
    if (pthread_mutex_unlock(&localtime_mutex) != 0)
    {
      r = NULL;
    }
  }
  return r;
}

static void
builtin_handler(SLEQP_LOG_LEVEL level, time_t time, const char* message)
{
  char buf[TIME_BUF_SIZE];

  struct tm result;

  localtime_synced(&time, &result);

  buf[strftime(buf, TIME_BUF_SIZE - 1, "%H:%M:%S", &result)] = '\0';

  fprintf(stderr,
          "[" SLEQP_FORMAT_BOLD "%s %s%5s" SLEQP_FORMAT_RESET "] %s\n",
          buf,
          level_infos[level].color,
          level_infos[level].name,
          message);
}

static SLEQP_LOG_HANDLER handler = builtin_handler;

void
sleqp_log_msg_level(int level, const char* fmt, ...)
{
  char message_buf[TOTAL_BUF_SIZE];

  time_t t = time(NULL);

  va_list args;

  va_start(args, fmt);
  vsnprintf(message_buf, TOTAL_BUF_SIZE, fmt, args);
  va_end(args);

  handler(level, t, message_buf);
}

void
sleqp_log_set_handler(SLEQP_LOG_HANDLER value)
{
  handler = value;
}

void
sleqp_log_trace_level(int level,
                      const char* file,
                      int line,
                      const char* fmt,
                      ...)
{
  char message_buf[TOTAL_BUF_SIZE];
  char total_buf[EXTENDED_BUF_SIZE];

  va_list args;

  va_start(args, fmt);
  vsnprintf(message_buf, TOTAL_BUF_SIZE, fmt, args);
  va_end(args);

  time_t t = time(NULL);

  snprintf(total_buf,
           EXTENDED_BUF_SIZE,
           SLEQP_FORMAT_DARK "%s:%d " SLEQP_FORMAT_RESET "%s",
           file,
           line,
           message_buf);

  handler(level, t, total_buf);
}
