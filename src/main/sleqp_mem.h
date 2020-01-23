#ifndef SLEQP_MEM_H
#define SLEQP_MEM_H

/**
 * @file sleqp_mem.h
 * @brief Definition of memory (de-)allocation functions.
 **/

#include <stdlib.h>

#include "sleqp_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define sleqp_allocate_memory(ptr, size)                                     \
  (*(ptr) = ((size) > 0) ? malloc((size)) : NULL), (((size) > 0) && (*(ptr) == NULL)) \
    ? SLEQP_NOMEM : SLEQP_OKAY

#define sleqp_reallocate_memory(ptr, size)                              \
  (*ptr = realloc(*ptr, size), *ptr != NULL) ? SLEQP_OKAY : SLEQP_NOMEM

  SLEQP_RETCODE sleqp_free(void** ptr);

#define sleqp_malloc(ptr)                       \
  sleqp_allocate_memory(ptr, sizeof(**ptr))

#define sleqp_calloc(ptr, count)                        \
  sleqp_allocate_memory(ptr, ((count) * sizeof(**ptr)))

#define sleqp_realloc(ptr, count)                               \
  sleqp_reallocate_memory(ptr, ((count) * sizeof(**ptr)))

#define sleqp_free(ptr)                         \
  free(*ptr); *ptr = NULL

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_MEM_H */
