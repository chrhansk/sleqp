#ifndef SLEQP_FEAS_H
#define SLEQP_FEAS_H

/**
 * @file sleqp_feas.h
 * @brief Definition of feasibility tests.
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

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_violated_constraint_multipliers(SleqpProblem* problem,
                                                      SleqpSparseVec* cons_vals,
                                                      SleqpSparseVec* multipliers,
                                                      SleqpWorkingSet* working_set);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_violated_variable_multipliers(SleqpProblem* problem,
                                                    SleqpSparseVec* primal,
                                                    SleqpSparseVec* multipliers,
                                                    SleqpWorkingSet* working_set);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_violated_constraints(SleqpProblem* problem,
                                           SleqpSparseVec* cons_val,
                                           int* violated_constraints,
                                           int* num_violated_constraints);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_violation_values(SleqpProblem* problem,
                                       SleqpIterate* iterate,
                                       SleqpSparseVec* violation);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_violation_inf_norm(SleqpProblem* problem,
                                         SleqpSparseVec* cons_val,
                                         double* max_violation);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_violation_one_norm(SleqpProblem* problem,
                                         SleqpSparseVec* cons_val,
                                         double* total_violation);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_FEAS_H */
