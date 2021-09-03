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

  SLEQP_RETCODE sleqp_dyn_func_val(SleqpFunc* func,
                                   double accuracy,
                                   double* func_val);

  SLEQP_RETCODE sleqp_dyn_func_cons_val(SleqpFunc* func,
                                        double accuracy,
                                        const SleqpSparseVec* cons_indices,
                                        SleqpSparseVec* cons_val);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_DYN_H */
