#ifndef SLEQP_TYPES_H
#define SLEQP_TYPES_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdbool.h>

#include "sleqp_log.h"

  enum SLEQP_Retcode
  {
    SLEQP_OKAY,
    SLEQP_NOMEM,
    SLEQP_INVALID,
    SLEQP_INTERNAL_ERROR,
  };

  typedef enum SLEQP_Retcode SLEQP_RETCODE;

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
      sleqp_log_error("Error <%d> in function call", status);              \
      return status;                                                       \
    }                                                                      \
  }                                                                        \
  while(0)

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_TYPES_H */
