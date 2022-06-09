#ifndef SLEQP_UTIL_H
#define SLEQP_UTIL_H

/**
 * @file util.h
 * @brief Definition of utility functions.
 **/

#include "func.h"
#include "iterate.h"
#include "problem.h"
#include "types.h"

#include "sparse/sparse_matrix.h"
#include "sparse/vec.h"

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_set_and_evaluate(SleqpProblem* problem,
                       SleqpIterate* iterate,
                       SLEQP_VALUE_REASON reason,
                       bool* reject);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_direction_in_working_set(SleqpProblem* problem,
                               const SleqpIterate* iterate,
                               const SleqpVec* direction,
                               double* cache,
                               double eps,
                               bool* in_working_set);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_max_step_length(const SleqpVec* x,
                      const SleqpVec* direction,
                      const SleqpVec* var_lb,
                      const SleqpVec* var_ub,
                      double* max_step_length);

double
sleqp_reduction_ratio(const double exact_reduction,
                      const double model_reduction);

#endif /* SLEQP_UTIL_H */
