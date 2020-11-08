#ifndef SLEQP_UTIL_H
#define SLEQP_UTIL_H

/**
 * @file sleqp_util.h
 * @brief Definition of utility functions.
 **/

#include "sleqp_func.h"
#include "sleqp_iterate.h"
#include "sleqp_problem.h"
#include "sleqp_types.h"

#include "sparse/sleqp_sparse_matrix.h"
#include "sparse/sleqp_sparse_vec.h"

#ifdef __cplusplus
extern "C" {
#endif

  SLEQP_RETCODE sleqp_set_and_evaluate(SleqpProblem* problem,
                                       SleqpIterate* iterate,
                                       SLEQP_VALUE_REASON reason);

  SLEQP_RETCODE sleqp_direction_in_working_set(SleqpProblem* problem,
                                               SleqpIterate* iterate,
                                               SleqpSparseVec* direction,
                                               double* cache,
                                               double eps,
                                               bool* in_working_set);

  SLEQP_RETCODE sleqp_max_step_length(SleqpSparseVec* x,
                                      SleqpSparseVec* direction,
                                      SleqpSparseVec* var_lb,
                                      SleqpSparseVec* var_ub,
                                      double* max_step_length);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_UTIL_H */
