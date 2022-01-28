#include "error.h"

#include <stdarg.h>
#include <threads.h>

#define ERROR_MSG_SIZE 2048

static thread_local char error_msg[ERROR_MSG_SIZE];
static thread_local SLEQP_ERROR_TYPE etype;

void
sleqp_error_type_set(SLEQP_ERROR_TYPE error_type)
{
  etype = error_type;
}

const char*
sleqp_error_msg()
{
  return error_msg;
}

SLEQP_ERROR_TYPE
sleqp_error_type()
{
  return etype;
}

void
sleqp_set_error(const char* file,
                int line,
                const char* func,
                SLEQP_ERROR_TYPE error_type,
                const char* fmt,
                ...)
{
  va_list args;

  etype = error_type;

  va_start(args, fmt);
  vsnprintf(error_msg, ERROR_MSG_SIZE, fmt, args);
  va_end(args);
}
