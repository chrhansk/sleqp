#ifndef SLEQP_PUB_ERROR_H
#define SLEQP_PUB_ERROR_H

#include "pub_types.h"

SLEQP_EXPORT SLEQP_ERROR_TYPE
sleqp_error_type();

SLEQP_EXPORT void
sleqp_error_type_set(SLEQP_ERROR_TYPE error_type);

SLEQP_EXPORT const char*
sleqp_error_msg();

SLEQP_EXPORT void
sleqp_set_error(const char* file,
                int line,
                const char* func,
                SLEQP_ERROR_TYPE error_type,
                const char* fmt,
                ...);

#define sleqp_raise(error_type, fmt, ...)                                      \
  do                                                                           \
  {                                                                            \
    sleqp_set_error(__FILE__,                                                  \
                    __LINE__,                                                  \
                    __PRETTY_FUNCTION__,                                       \
                    error_type,                                                \
                    fmt,                                                       \
                    ##__VA_ARGS__);                                            \
    return SLEQP_ERROR;                                                        \
  } while (false)

#endif /* SLEQP_PUB_ERROR_H */
