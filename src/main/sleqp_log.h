#ifndef SLEQP_LOG_H
#define SLEQP_LOG_H

#include <time.h>

#include "sleqp_defs.h"

/**
 * @file sleqp_log.h
 * @brief Definition of logging functions.
 **/

#if defined SLEQP_FORMAT_CODES

#define SLEQP_FORMAT_RESET "\x1B[0m"
#define SLEQP_FORMAT_RED "\x1B[31m"
#define SLEQP_FORMAT_GREEN "\x1B[32m"
#define SLEQP_FORMAT_YELLOW "\x1B[33m"
#define SLEQP_FORMAT_BLUE "\x1B[34m"
#define SLEQP_FORMAT_DARK "\x1b[90m"

#define SLEQP_FORMAT_BOLD "\x1B[1m"
#define SLEQP_FORMAT_NO_BOLD "\x1B[22m"

#else

#define SLEQP_FORMAT_RESET ""
#define SLEQP_FORMAT_RED ""
#define SLEQP_FORMAT_GREEN ""
#define SLEQP_FORMAT_YELLOW ""
#define SLEQP_FORMAT_BLUE ""
#define SLEQP_FORMAT_DARK ""

#define SLEQP_FORMAT_BOLD ""
#define SLEQP_FORMAT_NO_BOLD ""


#endif

#ifdef __cplusplus
extern "C" {
#endif

  typedef enum {
    SLEQP_LOG_ERROR = 0,
    SLEQP_LOG_WARN = 1,
    SLEQP_LOG_INFO = 2,
    SLEQP_LOG_DEBUG = 3,
    SLEQP_NUM_LOG_LEVELS = 4} SLEQP_LOG_LEVEL;

  SLEQP_LOG_LEVEL sleqp_log_level();

  void sleqp_log_set_level(SLEQP_LOG_LEVEL level);

  typedef void (*SLEQP_LOG_HANDLER)(SLEQP_LOG_LEVEL level,
                                    time_t time,
                                    const char* message);

  void sleqp_log_set_handler(SLEQP_LOG_HANDLER handler);

  void sleqp_log_msg_level(int level, const char *fmt, ...);

  void sleqp_log_trace_level(int level, const char *file, int line, const char *fmt, ...);

#define sleqp_log_log_trace(level, file, line, ...)           \
  do                                                          \
  {                                                           \
    if(sleqp_log_level() >= level)                            \
    {                                                         \
      sleqp_log_trace_level(level, file, line, __VA_ARGS__);  \
    }                                                         \
  } while(0)

#define sleqp_log_log_msg(level, ...)           \
  do                                            \
  {                                             \
    if(sleqp_log_level() >= level)              \
    {                                           \
      sleqp_log_msg_level(level, __VA_ARGS__);  \
    }                                           \
  } while(0)

#define sleqp_log_info(...)  sleqp_log_log_msg(SLEQP_LOG_INFO, __VA_ARGS__)
#define sleqp_log_warn(...)  sleqp_log_log_msg(SLEQP_LOG_WARN, __VA_ARGS__)
#define sleqp_log_error(...) sleqp_log_log_msg(SLEQP_LOG_ERROR, __VA_ARGS__)

#define sleqp_log_debug(...) sleqp_log_log_msg(SLEQP_LOG_DEBUG, __VA_ARGS__)

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_LOG_H */
