#ifndef SLEQP_DYN_H
#define SLEQP_DYN_H

#include "pub_dyn.h"

SLEQP_RETCODE
sleqp_dyn_func_set_accuracy(SleqpFunc* func, double accuracy);

SLEQP_RETCODE
sleqp_dyn_func_get_accuracy(SleqpFunc* func, double* accuracy);

SLEQP_RETCODE
sleqp_dyn_func_val(SleqpFunc* func, double accuracy, double* func_val);

SLEQP_RETCODE
sleqp_dyn_func_cons_val(SleqpFunc* func,
                        double accuracy,
                        const SleqpSparseVec* cons_indices,
                        SleqpSparseVec* cons_val);

#endif /* SLEQP_DYN_H */
