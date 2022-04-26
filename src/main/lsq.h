#ifndef SLEQP_LSQ_H
#define SLEQP_LSQ_H

#include "pub_lsq.h"

#include "func.h"
#include "timer.h"

double
sleqp_lsq_func_get_levenberg_marquardt(SleqpFunc* func);

int
sleqp_lsq_func_num_residuals(SleqpFunc* func);

SLEQP_NODISCARD SLEQP_RETCODE
sleqp_lsq_func_nonzeros(SleqpFunc* func,
                        int* residual_nnz,
                        int* jac_fwd_nnz,
                        int* jac_adj_nnz,
                        int* cons_val_nnz,
                        int* cons_jac_nnz);

SLEQP_NODISCARD SLEQP_RETCODE
sleqp_lsq_func_residuals(SleqpFunc* func, SleqpVec* residuals);

SLEQP_NODISCARD SLEQP_RETCODE
sleqp_lsq_func_jac_forward(SleqpFunc* func,
                           const SleqpVec* forward_direction,
                           SleqpVec* product);

SLEQP_NODISCARD SLEQP_RETCODE
sleqp_lsq_func_jac_adjoint(SleqpFunc* func,
                           const SleqpVec* adjoint_direction,
                           SleqpVec* product);

SleqpTimer*
sleqp_lsq_func_residual_timer(SleqpFunc* func);

SleqpTimer*
sleqp_lsq_func_adjoint_timer(SleqpFunc* func);

SleqpTimer*
sleqp_lsq_func_forward_timer(SleqpFunc* func);

void*
sleqp_lsq_func_get_data(SleqpFunc* func);

#endif /* SLEQP_LSQ_H */
