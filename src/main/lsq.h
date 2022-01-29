#ifndef SLEQP_LSQ_H
#define SLEQP_LSQ_H

#include "pub_lsq.h"

#include "func.h"

double
sleqp_lsq_func_get_levenberg_marquardt(SleqpFunc* func);

int
sleqp_lsq_func_num_residuals(SleqpFunc* func);

SLEQP_NODISCARD SLEQP_RETCODE
sleqp_lsq_func_residuals(SleqpFunc* func, SleqpSparseVec* residuals);

SLEQP_NODISCARD SLEQP_RETCODE
sleqp_lsq_func_jac_forward(SleqpFunc* func,
                           const SleqpSparseVec* forward_direction,
                           SleqpSparseVec* product);

SLEQP_NODISCARD SLEQP_RETCODE
sleqp_lsq_func_jac_adjoint(SleqpFunc* func,
                           const SleqpSparseVec* adjoint_direction,
                           SleqpSparseVec* product);

void*
sleqp_lsq_func_get_data(SleqpFunc* func);

#endif /* SLEQP_LSQ_H */
