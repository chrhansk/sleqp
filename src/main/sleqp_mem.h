#ifndef SLEQP_MEM_H
#define SLEQP_MEM_H

#include <stdlib.h>

#include "sleqp_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define sleqp_allocate_memory(ptr, size)        \
  {                                             \
    *ptr = malloc(size);                        \
                                                \
    if(!*ptr)                                   \
    {                                           \
      return SLEQP_NOMEM;                       \
    }                                           \
  }

  SLEQP_RETCODE sleqp_free(void** ptr);

#define sleqp_malloc(ptr)                       \
  do                                            \
  {                                             \
    sleqp_allocate_memory(ptr,                  \
                          sizeof(**ptr));       \
  }                                             \
  while(0)

#define sleqp_calloc(ptr, count)                        \
  do                                                    \
  {                                                     \
    sleqp_allocate_memory(ptr,                          \
                          (count)*sizeof(**ptr));       \
  }                                                     \
  while(0)

#define sleqp_realloc(ptr, count)                       \
  do                                                    \
  {                                                     \
    *ptr = realloc(ptr,                                 \
                  (count)*sizeof(**ptr));               \
    if(!*ptr)                                           \
    {                                                   \
      return SLEQP_NOMEM;                               \
    }                                                   \
  }                                                     \
  while(0)

#define sleqp_free(ptr)                         \
  do                                            \
  {                                             \
    free(ptr);                                  \
  }                                             \
  while(0)

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_MEM_H */
