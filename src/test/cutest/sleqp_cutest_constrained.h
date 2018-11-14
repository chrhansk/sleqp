#ifndef SLEQP_CUTEST_CONSTRAINED_H
#define SLEQP_CUTEST_CONSTRAINED_H

#include "sleqp.h"

#ifdef __cplusplus
extern "C" {
#endif

  SLEQP_RETCODE sleqp_cutest_cons_func_create(SleqpFunc** star,
                                              int num_variables,
                                              int num_constraints,
                                              double eps);

  SLEQP_RETCODE sleqp_cutest_cons_func_free(SleqpFunc** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_CUTEST_CONSTRAINED_H */
