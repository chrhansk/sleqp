#include "sleqp_callback_handler.h"

#include <assert.h>

#include "sleqp_cmp.h"
#include "sleqp_mem.h"

typedef struct Callback
{
  void* callback;
  void* callback_data;
} Callback;

struct SleqpCallbackHandler
{
  int refcount;

  int size;
  int capacity;

  Callback* callbacks;
};

SLEQP_RETCODE sleqp_callback_handler_create(SleqpCallbackHandler** star)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpCallbackHandler* handler = *star;

  *handler = (SleqpCallbackHandler){0};

  handler->refcount = 1;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_callback_handler_capture(SleqpCallbackHandler* handler)
{
  ++handler->refcount;

  return SLEQP_OKAY;
}

int sleqp_callback_handler_size(SleqpCallbackHandler* handler)
{
  return handler->size;
}

SLEQP_RETCODE sleqp_callback_handler_get(SleqpCallbackHandler* handler,
                                         int pos,
                                         void** callback,
                                         void** callback_data)
{
  if(pos < 0 || pos >= handler->size)
  {
    sleqp_log_error("Invalid index");
    return SLEQP_ILLEGAL_ARGUMENT;
  }

  assert(callback);
  assert(callback_data);

  (*callback) = handler->callbacks[pos].callback;
  (*callback_data) = handler->callbacks[pos].callback_data;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE callback_handler_reserve(SleqpCallbackHandler* handler)
{
  if(handler->size < handler->capacity)
  {
    return SLEQP_OKAY;
  }

  const int next_capacity = SLEQP_MAX(1, 2*handler->capacity);

  SLEQP_CALL(sleqp_realloc(&(handler->callbacks), next_capacity));

  handler->capacity = next_capacity;

  assert(handler->size < handler->capacity);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_callback_handler_add(SleqpCallbackHandler* handler,
                                         void* callback,
                                         void* callback_data)
{
  SLEQP_CALL(callback_handler_reserve(handler));

  handler->callbacks[handler->size] = (Callback) {
    .callback = callback,
    .callback_data = callback_data
  };

  ++(handler->size);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_callback_handler_remove(SleqpCallbackHandler* handler,
                                            void* callback,
                                            void* callback_data)
{
  bool found = false;

  int pos;

  for(pos = 0; pos < handler->size;++pos)
  {
    if((handler->callbacks[pos].callback == callback) &&
       (handler->callbacks[pos].callback_data == callback_data))
    {
      found = true;
      break;
    }
  }

  if(!found)
  {
    sleqp_log_error("Attempted to remove non-existing callback");
    return SLEQP_ILLEGAL_ARGUMENT;
  }

  --handler->size;

  // not last element
  if(pos != handler->size)
  {
    handler->callbacks[pos] = handler->callbacks[handler->size];
  }

  handler->callbacks[handler->size] = (Callback) {0};

  return SLEQP_OKAY;
}

static SLEQP_RETCODE callback_handler_free(SleqpCallbackHandler** star)
{
  SleqpCallbackHandler* handler = *star;

  if(!handler)
  {
    return SLEQP_OKAY;
  }

  sleqp_free(&(handler->callbacks));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_callback_handler_release(SleqpCallbackHandler** star)
{
  SleqpCallbackHandler* handler = *star;

  if(!handler)
  {
    return SLEQP_OKAY;
  }

  if(--handler->refcount == 0)
  {
    SLEQP_CALL(callback_handler_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
