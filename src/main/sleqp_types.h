#ifndef SLEQP_TYPES_H
#define SLEQP_TYPES_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>

  enum SLEQP_Retcode
  {
    SLEQP_OKAY,
    SLEQP_NOMEM,
    SLEQP_INVALID,
  };

  typedef enum SLEQP_Retcode SLEQP_RETCODE;

  typedef unsigned int SLEQP_Bool;

#define SLEQP_CALL(x)                                                   \
  do                                                                    \
  {                                                                     \
    SLEQP_RETCODE _retcode_;                                            \
    if( (_retcode_ = (x)) != SLEQP_OKAY )                               \
    {                                                                   \
      fprintf(stderr, "Error <%d> in function call\n", _retcode_);      \
      return _retcode_;                                                 \
    }                                                                   \
  }                                                                     \
  while(0)

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_TYPES_H */
