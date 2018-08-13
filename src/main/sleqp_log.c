#include "sleqp_log.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>

#define TIME_BUF_SIZE 128

#define NORM  "\x1B[0m"
#define RED  "\x1B[31m"
#define GREEN  "\x1B[32m"
#define YELLOW  "\x1B[33m"
#define BLUE  "\x1B[34m"
#define DARK "\x1b[90m"

#define BOLD "\x1B[1m"
#define NO_BOLD "\x1B[22m"

struct LevelInfo
{
  const char* name;
  const char* color;
};

static struct LevelInfo level_infos[SLEQP_NUM_LOG_LEVELS] =
{
  {"error", RED},
  {"warn", YELLOW},
  {"info", GREEN},
  {"debug", BLUE}
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
          "[" BOLD "%s %s%s" NORM NO_BOLD "] " NORM,
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
          "[" BOLD "%s %s%s" NORM NO_BOLD "] " DARK "%s:%d " NORM,
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
