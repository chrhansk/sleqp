#ifndef SLEQP_LOG_H
#define SLEQP_LOG_H

#define SLEQP_NORM  "\x1B[0m"
#define SLEQP_RED  "\x1B[31m"
#define SLEQP_GREEN  "\x1B[32m"
#define SLEQP_YELLOW  "\x1B[33m"
#define SLEQP_BLUE  "\x1B[34m"
#define SLEQP_DARK "\x1b[90m"

#define SLEQP_BOLD "\x1B[1m"
#define SLEQP_NO_BOLD "\x1B[22m"

#ifdef __cplusplus
extern "C" {
#endif

enum SLEQP_Log_Level {SLEQP_LOG_ERROR = 0,
                      SLEQP_LOG_WARN = 1,
                      SLEQP_LOG_INFO = 2,
                      SLEQP_LOG_DEBUG = 3,
                      SLEQP_NUM_LOG_LEVELS = 4};

typedef enum SLEQP_Log_Level SLEQP_LOG_LEVEL;

SLEQP_LOG_LEVEL sleqp_log_level();

void sleqp_log_msg_level(int level, const char *fmt, ...);

void sleqp_log_trace_level(int level, const char *file, int line, const char *fmt, ...);

#define sleqp_log_log_trace(level, file, line, ...)                     \
  do                                                                    \
  {                                                                     \
    if(sleqp_log_level() >= level)                                      \
  {                                                                     \
    sleqp_log_trace_level(level, file, line, __VA_ARGS__);              \
  }                                                                     \
  } while(0)

#define sleqp_log_log_msg(level, ...)                           \
  do                                                            \
  {                                                             \
    if(sleqp_log_level() >= level)                              \
    {                                                           \
      sleqp_log_msg_level(level, __VA_ARGS__);                  \
    }                                                           \
  } while(0)

#define sleqp_log_info(...)  sleqp_log_log_msg(SLEQP_LOG_INFO, __VA_ARGS__)

#define sleqp_log_debug(...) sleqp_log_log_trace(SLEQP_LOG_DEBUG, __FILE__, __LINE__, __VA_ARGS__)
#define sleqp_log_warn(...)  sleqp_log_log_trace(SLEQP_LOG_WARN,  __FILE__, __LINE__, __VA_ARGS__)
#define sleqp_log_error(...) sleqp_log_log_trace(SLEQP_LOG_ERROR,  __FILE__, __LINE__, __VA_ARGS__)

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_LOG_H */
