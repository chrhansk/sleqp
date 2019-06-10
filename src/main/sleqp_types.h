#ifndef SLEQP_TYPES_H
#define SLEQP_TYPES_H

/**
 * @file sleqp_types.h
 * @brief Definition of basic types.
 **/

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdbool.h>

#include "sleqp_log.h"

#define RETCODES C(SLEQP_OKAY, 0)               \
    C(SLEQP_NOMEM, 1)                           \
    C(SLEQP_ILLEGAL_ARGUMENT, 2)                \
    C(SLEQP_INVALID_DERIV, 3)                   \
    C(SLEQP_INTERNAL_ERROR, 4)                  \
    C(SLEQP_NUM_RETCODES, 5)

#define C(k, v) k = v,
  typedef enum { RETCODES } SLEQP_RETCODE;
#undef C

#define C(k, v) [v] = #k,
  static const char * const sleqp_retcode_names[] = { RETCODES };
#undef C
#undef RETCODES

  typedef enum {SLEQP_INACTIVE = 0,
                SLEQP_ACTIVE_LOWER = (1 << 1),
                SLEQP_ACTIVE_UPPER = (1 << 2),
                SLEQP_ACTIVE_BOTH  = (1 << 3)} SLEQP_ACTIVE_STATE;

  typedef enum {SLEQP_OPTIMAL,
                SLEQP_FEASIBLE,
                SLEQP_INFEASIBLE,
                SLEQP_INVALID} SLEQP_STATUS;

#define SLEQP_CALL(x)                                                      \
  do                                                                       \
  {                                                                        \
    SLEQP_RETCODE status;                                                  \
    if( (status = (x)) != SLEQP_OKAY )                                     \
    {                                                                      \
       sleqp_log_error("Error <%d> (%s) in function %s",                   \
                       status,                                             \
                       sleqp_retcode_names[status],                        \
                       __func__);                                          \
      return status;                                                       \
    }                                                                      \
  }                                                                        \
  while(0)

  typedef enum {
    SLEQP_STEPTYPE_NONE = 0,
    SLEQP_STEPTYPE_ACCEPTED,
    SLEQP_STEPTYPE_ACCEPTED_FULL,
    SLEQP_STEPTYPE_SOC_ACCEPTED,
    SLEQP_STEPTYPE_REJECTED
  } SLEQP_STEPTYPE;

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_TYPES_H */
