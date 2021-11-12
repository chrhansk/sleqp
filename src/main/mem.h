#ifndef SLEQP_MEM_H
#define SLEQP_MEM_H

/**
 * @file mem.h
 * @brief Definition of memory (de-)allocation functions.
 **/

#include <stdlib.h>

#include "types.h"

#define sleqp_allocate_memory(ptr, size)                                \
        (*(ptr) = ((size) > 0) ? malloc((size)) : NULL), (((size) > 0) && (*(ptr) == NULL)) \
                ? SLEQP_NOMEM : SLEQP_OKAY

#define sleqp_reallocate_memory(ptr, size)                              \
        (*ptr = realloc(*ptr, size), (((size) > 0) && (*(ptr) == NULL))) ? SLEQP_NOMEM : SLEQP_OKAY

SLEQP_NODISCARD SLEQP_RETCODE sleqp_free(void** ptr);

#define sleqp_malloc(ptr)                               \
        sleqp_allocate_memory(ptr, sizeof(**ptr))

#define sleqp_alloc_array(ptr, count)                           \
        sleqp_allocate_memory(ptr, ((count) * sizeof(**ptr)))

#define sleqp_realloc(ptr, count)                               \
        sleqp_reallocate_memory(ptr, ((count) * sizeof(**ptr)))

#define sleqp_free(ptr)                         \
        free(*ptr); *ptr = NULL

#endif /* SLEQP_MEM_H */
