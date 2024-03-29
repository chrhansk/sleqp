#ifndef SLEQP_CALLBACK_HANDLER_H
#define SLEQP_CALLBACK_HANDLER_H

#include "types.h"

typedef struct SleqpCallbackHandler SleqpCallbackHandler;

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_callback_handler_create(SleqpCallbackHandler** star);

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_callback_handler_capture(SleqpCallbackHandler* handler);

int
sleqp_callback_handler_size(SleqpCallbackHandler* handler);

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_callback_handler_get(SleqpCallbackHandler* handler,
                           int pos,
                           void** callback,
                           void** callback_data);

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_callback_handler_add(SleqpCallbackHandler* handler,
                           void* callback,
                           void* callback_data);

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_callback_handler_remove(SleqpCallbackHandler* handler,
                              void* callback,
                              void* callback_data);

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_callback_handler_release(SleqpCallbackHandler** star);

#define SLEQP_CALLBACK_HANDLER_EXECUTE(handler, func_type, ...)                \
  do                                                                           \
  {                                                                            \
    const int size = sleqp_callback_handler_size(handler);                     \
                                                                               \
    void* callback_func;                                                       \
    void* callback_data;                                                       \
                                                                               \
    for (int pos = 0; pos < size; ++pos)                                       \
    {                                                                          \
      SLEQP_CALL(sleqp_callback_handler_get(handler,                           \
                                            pos,                               \
                                            &callback_func,                    \
                                            &callback_data));                  \
                                                                               \
      SLEQP_CALL(((func_type)(callback_func))(__VA_ARGS__, callback_data));    \
    }                                                                          \
  } while (false)

#endif /* SLEQP_CALLBACK_HANDLER_H */
