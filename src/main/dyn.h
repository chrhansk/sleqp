#ifndef SLEQP_DYN_H
#define SLEQP_DYN_H

#include "pub_dyn.h"

SLEQP_WARNUNUSED SLEQP_RETCODE
sleqp_dyn_func_eval(SleqpFunc* func,
                    double* obj_val,
                    SleqpVec* cons_val,
                    double* error);

SLEQP_WARNUNUSED SLEQP_RETCODE
sleqp_dyn_func_set_error_bound(SleqpFunc* func, double error_bound);

SLEQP_WARNUNUSED SLEQP_RETCODE
sleqp_dyn_func_error(SleqpFunc* func, double* error);

SLEQP_WARNUNUSED SLEQP_RETCODE
sleqp_dyn_func_error_bound(SleqpFunc* func, double* error_bound);

SLEQP_WARNUNUSED SLEQP_RETCODE
sleqp_dyn_func_error_estimate(SleqpFunc* func, double* error_estimate);

SLEQP_WARNUNUSED SLEQP_RETCODE
sleqp_dyn_func_set_obj_weight(SleqpFunc* func, double obj_weight);

SLEQP_WARNUNUSED SLEQP_RETCODE
sleqp_dyn_func_set_cons_weights(SleqpFunc* func, const double* cons_weights);

SLEQP_WARNUNUSED SLEQP_RETCODE
sleqp_dyn_set_penalty_cons_weights(SleqpFunc* func,
                                   double penalty,
                                   double* cache);

#endif /* SLEQP_DYN_H */
