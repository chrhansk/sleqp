#ifndef SLEQP_AMPL_MEM_H
#define SLEQP_AMPL_MEM_H

/**
 * @file mem.h
 * @brief Definition of memory (de-)allocation functions.
 **/

#include <stdlib.h>

#include "sleqp/pub_types.h"

#ifdef __cplusplus
extern "C"
{
#endif

#define sleqp_ampl_allocate_memory(ptr, size)                                         \
  (*(ptr) = ((size) > 0) ? malloc((size)) : NULL), (((size) > 0) && (*(ptr) == NULL)) \
                                                       ? SLEQP_NOMEM                  \
                                                       : SLEQP_OKAY

#define sleqp_ampl_reallocate_memory(ptr, size) \
  (*ptr = realloc(*ptr, size), (((size) > 0) && (*(ptr) == NULL))) ? SLEQP_NOMEM : SLEQP_OKAY

  SLEQP_NODISCARD SLEQP_RETCODE sleqp_ampl_free(void **ptr);

#define sleqp_ampl_malloc(ptr) \
  sleqp_ampl_allocate_memory(ptr, sizeof(**ptr))

#define sleqp_ampl_alloc_array(ptr, count) \
  sleqp_ampl_allocate_memory(ptr, ((count) * sizeof(**ptr)))

#define sleqp_ampl_realloc(ptr, count) \
  sleqp_ampl_reallocate_memory(ptr, ((count) * sizeof(**ptr)))

#define sleqp_ampl_free(ptr) \
  free(*ptr);                \
  *ptr = NULL

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_AMPL_MEM_H */
