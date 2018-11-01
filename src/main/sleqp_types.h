#ifndef SLEQP_TYPES_H
#define SLEQP_TYPES_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdbool.h>

#include "sleqp_log.h"

#define RETCODES C(SLEQP_OKAY, 0)               \
    C(SLEQP_NOMEM, 1)                           \
    C(SLEQP_INVALID, 2)                         \
    C(SLEQP_INVALID_DERIV, 3)                   \
    C(SLEQP_INTERNAL_ERROR, 4)                  \
    C(SLEQP_NUM_RETCODES, 5)

#define C(k, v) k = v,
  enum SLEQP_Retcode { RETCODES };
#undef C

  typedef enum SLEQP_Retcode SLEQP_RETCODE;

#define C(k, v) [v] = #k,
  static const char * const sleqp_retcode_names[] = { RETCODES };
#undef C
#undef RETCODES

  enum SLEQP_Active_State {SLEQP_INACTIVE = 0,
                           SLEQP_ACTIVE_LOWER = (1 << 1),
                           SLEQP_ACTIVE_UPPER = (1 << 2),
                           SLEQP_ACTIVE = SLEQP_ACTIVE_LOWER | SLEQP_ACTIVE_UPPER};

  typedef enum SLEQP_Active_State SLEQP_ACTIVE_STATE;

#define SLEQP_CALL(x)                                                      \
  do                                                                       \
  {                                                                        \
    SLEQP_RETCODE status;                                                  \
    if( (status = (x)) != SLEQP_OKAY )                                     \
    {                                                                      \
      sleqp_log_error("Error <%d> (%s) in function call",                  \
                      status,                                              \
                      sleqp_retcode_names[status]);                        \
      return status;                                                       \
    }                                                                      \
  }                                                                        \
  while(0)

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_TYPES_H */
