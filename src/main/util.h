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
#include "sparse/sparse_vec.h"

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_set_and_evaluate(SleqpProblem* problem,
                                     SleqpIterate* iterate,
                                     SLEQP_VALUE_REASON reason,
                                     bool* reject);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_direction_in_working_set(SleqpProblem* problem,
                                             const SleqpIterate* iterate,
                                             const SleqpSparseVec* direction,
                                             double* cache,
                                             double eps,
                                             bool* in_working_set);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_max_step_length(const SleqpSparseVec* x,
                                    const SleqpSparseVec* direction,
                                    const SleqpSparseVec* var_lb,
                                    const SleqpSparseVec* var_ub,
                                    double* max_step_length);

#endif /* SLEQP_UTIL_H */
