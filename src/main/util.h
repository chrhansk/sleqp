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

#ifdef __cplusplus
extern "C" {
#endif

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_set_and_evaluate(SleqpProblem* problem,
                                       SleqpIterate* iterate,
                                       SLEQP_VALUE_REASON reason);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_direction_in_working_set(SleqpProblem* problem,
                                               SleqpIterate* iterate,
                                               SleqpSparseVec* direction,
                                               double* cache,
                                               double eps,
                                               bool* in_working_set);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_max_step_length(SleqpSparseVec* x,
                                      SleqpSparseVec* direction,
                                      SleqpSparseVec* var_lb,
                                      SleqpSparseVec* var_ub,
                                      double* max_step_length);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_UTIL_H */
