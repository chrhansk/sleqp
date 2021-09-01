#ifndef SLEQP_DYN_H
#define SLEQP_DYN_H

#include "pub_dyn.h"

#ifdef __cplusplus
extern "C" {
#endif

  SLEQP_RETCODE sleqp_dyn_func_set_accuracy(SleqpFunc* func,
                                            double accuracy);

  SLEQP_RETCODE sleqp_dyn_func_get_accuracy(SleqpFunc* func,
                                            double* accuracy);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_DYN_H */
