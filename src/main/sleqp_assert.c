#include "sleqp_assert.h"

#include <assert.h>
#include <stdarg.h>

#include "sleqp_log.h"

#define BUF_SIZE 2048

void sleqp_log_assert_fail(const char *assertion,
                           const char *file,
                           unsigned int line,
                           const char *function)
{
  sleqp_log_error("%s:%d: %s: Assertion `%s' failed",
                  file,
                  line,
                  function,
                  assertion);

  assert(false);
}

void sleqp_log_assert_fail_msg(const char *file,
                               unsigned int line,
                               const char *function,
                               const char* format,
                               ...)
{
  char message_buf[BUF_SIZE];

  va_list args;

  va_start(args, format);
  vsnprintf(message_buf, BUF_SIZE, format, args);
  va_end(args);

  sleqp_log_error("%s:%d: %s: %s",
                  file,
                  line,
                  function,
                  message_buf);

  assert(false);
}
