#ifndef SLEQP_FEAS_H
#define SLEQP_FEAS_H

/**
 * @file feas.h
 * @brief Definition of feasibility tests.
 **/

#include "func.h"
#include "iterate.h"
#include "problem.h"
#include "types.h"

#include "sparse/sparse_matrix.h"
#include "sparse/sparse_vec.h"

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_violated_constraint_multipliers(SleqpProblem* problem,
                                                    const SleqpSparseVec* cons_vals,
                                                    SleqpSparseVec* multipliers,
                                                    SleqpWorkingSet* working_set);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_violated_constraints(SleqpProblem* problem,
                                         SleqpSparseVec* cons_val,
                                         int* violated_constraints,
                                         int* num_violated_constraints);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_violation_values(SleqpProblem* problem,
                                     const SleqpSparseVec* cons_val,
                                     SleqpSparseVec* violation);

/**
 * Computes the residuals of the given constraint values with respect
 * to the upper / lower bounds of the underlying problem.
 * The residuals are unsigned (i.e., always non-negative).
 **/
SLEQP_NODISCARD
SLEQP_RETCODE sleqp_feasibility_residuals(SleqpProblem* problem,
                                          const SleqpSparseVec* cons_val,
                                          SleqpSparseVec* residuals,
                                          SleqpWorkingSet* working_set);

/**
 * Computes the residuals of the given constraint values with respect
 * to the upper / lower bounds of the underlying problem.
 * The residuals are signed, i.e., positive if a constraint value exceeds its
 * upper bound, and negative if the lower bound exceeds the constraint value.
 **/
SLEQP_NODISCARD
SLEQP_RETCODE sleqp_signed_feasibility_residuals(SleqpProblem* problem,
                                                 const SleqpSparseVec* cons_val,
                                                 SleqpSparseVec* residuals,
                                                 SleqpWorkingSet* working_set);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_violation_inf_norm(SleqpProblem* problem,
                                       SleqpSparseVec* cons_val,
                                       double* max_violation);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_violation_one_norm(SleqpProblem* problem,
                                       SleqpSparseVec* cons_val,
                                       double* total_violation);

#endif /* SLEQP_FEAS_H */
