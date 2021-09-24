#ifndef ZERO_FUNC_H
#define ZERO_FUNC_H

#include "func.h"

#ifdef __cplusplus
extern "C" {
#endif

  SLEQP_RETCODE zero_func_create(SleqpFunc** star,
                                 int num_variables,
                                 int num_constraints);

#ifdef __cplusplus
}
#endif

#endif /* ZERO_FUNC_H */
