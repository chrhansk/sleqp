#ifndef SLEQP_PUB_LOG_H
#define SLEQP_PUB_LOG_H

#include <time.h>

#include "defs.h"
#include "export.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef enum {
    SLEQP_LOG_ERROR = 0,
    SLEQP_LOG_WARN = 1,
    SLEQP_LOG_INFO = 2,
    SLEQP_LOG_DEBUG = 3,
    SLEQP_NUM_LOG_LEVELS = 4}
    SLEQP_LOG_LEVEL;

  SLEQP_EXPORT SLEQP_LOG_LEVEL sleqp_log_level();

  SLEQP_EXPORT void sleqp_log_set_level(SLEQP_LOG_LEVEL level);

  typedef void (*SLEQP_LOG_HANDLER)(SLEQP_LOG_LEVEL level,
                                    time_t time,
                                    const char* message);

  SLEQP_EXPORT void sleqp_log_set_handler(SLEQP_LOG_HANDLER handler);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_PUB_LOG_H */
