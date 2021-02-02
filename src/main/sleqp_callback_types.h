#ifndef SLEQP_CALLBACK_H
#define SLEQP_CALLBACK_H

#include "sleqp_types.h"
#include "sleqp_iterate.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpSolver SleqpSolver;

  typedef SLEQP_RETCODE (*SLEQP_ACCEPTED_ITERATE)(SleqpSolver* solver,
                                                  SleqpIterate* iterate,
                                                  void* callback_data);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_CALLBACK_H */
