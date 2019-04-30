#include "sleqp_log.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>

#define TIME_BUF_SIZE 128

struct LevelInfo
{
  const char* name;
  const char* color;
};

static struct LevelInfo level_infos[SLEQP_NUM_LOG_LEVELS] =
{
  {"error", SLEQP_FORMAT_RED},
  {"warn", SLEQP_FORMAT_YELLOW},
  {"info", SLEQP_FORMAT_GREEN},
  {"debug", SLEQP_FORMAT_BLUE}
};

static SLEQP_LOG_LEVEL level = SLEQP_LOG_DEBUG;

SLEQP_LOG_LEVEL sleqp_log_level()
{
  return level;
}

void sleqp_log_msg_level(int level, const char *fmt, ...)
{

  char buf[TIME_BUF_SIZE];

  time_t t = time(NULL);
  struct tm *lt = localtime(&t);

  buf[strftime(buf, TIME_BUF_SIZE - 1, "%H:%M:%S", lt)] = '\0';

  fprintf(stderr,
          "[" SLEQP_FORMAT_BOLD "%s %s%5s" SLEQP_FORMAT_RESET "] ",
          buf,
          level_infos[level].color,
          level_infos[level].name);

  va_list args;

  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);

  fprintf(stderr, "\n");

}


void sleqp_log_trace_level(int level, const char *file, int line, const char *fmt, ...)
{

  char buf[TIME_BUF_SIZE];

  time_t t = time(NULL);
  struct tm *lt = localtime(&t);

  buf[strftime(buf, TIME_BUF_SIZE - 1, "%H:%M:%S", lt)] = '\0';

  fprintf(stderr,
          "[" SLEQP_FORMAT_BOLD "%s %s%5s" SLEQP_FORMAT_RESET "] "
          SLEQP_FORMAT_DARK "%s:%d " SLEQP_FORMAT_RESET,
          buf,
          level_infos[level].color,
          level_infos[level].name,
          file,
          line);

  va_list args;

  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);

  fprintf(stderr, "\n");

}
