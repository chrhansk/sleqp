#ifndef SLEQP_CUTEST_CONSTRAINED_H
#define SLEQP_CUTEST_CONSTRAINED_H

#include "sleqp.h"

#include "sleqp_cutest_data.h"
#include "sleqp_cutest_types.h"

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_cutest_cons_func_create(SleqpFunc** star,
                              int num_variables,
                              int num_constraints,
                              int num_linear,
                              SleqpParams* params);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_cutest_cons_problem_create(SleqpProblem** star,
                                 SleqpCutestData* data,
                                 SleqpParams* params,
                                 bool force_nonlinear);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_cutest_eval_linear(SleqpFunc* func, SleqpVec* linear);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_cutest_eval_linear_coeffs(SleqpFunc* func, SleqpMat* coeffs);

#endif /* SLEQP_CUTEST_CONSTRAINED_H */
