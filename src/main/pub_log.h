#ifndef SLEQP_PUB_LOG_H
#define SLEQP_PUB_LOG_H

#include <time.h>

#include "sleqp/defs.h"
#include "sleqp/export.h"

typedef enum
{
  SLEQP_LOG_SILENT     = 0,
  SLEQP_LOG_ERROR      = 1,
  SLEQP_LOG_WARN       = 2,
  SLEQP_LOG_INFO       = 3,
  SLEQP_LOG_DEBUG      = 4,
  SLEQP_NUM_LOG_LEVELS = 5
} SLEQP_LOG_LEVEL;

SLEQP_EXPORT SLEQP_LOG_LEVEL
sleqp_log_level();

SLEQP_EXPORT void
sleqp_log_set_level(SLEQP_LOG_LEVEL level);

typedef void (*SLEQP_LOG_HANDLER)(SLEQP_LOG_LEVEL level,
                                  time_t time,
                                  const char* message);

SLEQP_EXPORT void
sleqp_log_set_handler(SLEQP_LOG_HANDLER handler);

SLEQP_EXPORT void
sleqp_log_msg_level(int level, const char* fmt, ...) SLEQP_FORMAT_PRINTF(2, 3);

void
sleqp_log_trace_level(int level,
                      const char* file,
                      int line,
                      const char* fmt,
                      ...) SLEQP_FORMAT_PRINTF(4, 5);

#define sleqp_bool_string(x) ((x) ? "true" : "false")

#define sleqp_log_log_trace(level, file, line, ...)                            \
  do                                                                           \
  {                                                                            \
    if (sleqp_log_level() >= level)                                            \
    {                                                                          \
      sleqp_log_trace_level(level, file, line, __VA_ARGS__);                   \
    }                                                                          \
  } while (0)

#define sleqp_log_log_msg(level, ...)                                          \
  do                                                                           \
  {                                                                            \
    if (sleqp_log_level() >= level)                                            \
    {                                                                          \
      sleqp_log_msg_level(level, __VA_ARGS__);                                 \
    }                                                                          \
  } while (0)

#define sleqp_log_info(...) sleqp_log_log_msg(SLEQP_LOG_INFO, __VA_ARGS__)
#define sleqp_log_warn(...) sleqp_log_log_msg(SLEQP_LOG_WARN, __VA_ARGS__)
#define sleqp_log_error(...) sleqp_log_log_msg(SLEQP_LOG_ERROR, __VA_ARGS__)

#define sleqp_log_debug(...) sleqp_log_log_msg(SLEQP_LOG_DEBUG, __VA_ARGS__)

#endif /* SLEQP_PUB_LOG_H */
