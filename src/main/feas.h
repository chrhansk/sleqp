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

#include "sparse/mat.h"
#include "sparse/vec.h"

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_violated_cons_multipliers(SleqpProblem* problem,
                                const SleqpVec* cons_vals,
                                SleqpVec* multipliers,
                                SleqpWorkingSet* working_set);

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_violated_constraints(SleqpProblem* problem,
                           SleqpVec* cons_val,
                           int* violated_constraints,
                           int* num_violated_constraints);

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_violation_values(SleqpProblem* problem,
                       const SleqpVec* cons_val,
                       SleqpVec* violation);

/**
 * Computes the residuals of the given constraint values with respect
 * to the upper / lower bounds of the underlying problem.
 * The residuals are unsigned (i.e., always non-negative).
 **/
SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_feasibility_residuals(SleqpProblem* problem,
                            const SleqpVec* cons_val,
                            SleqpVec* residuals,
                            SleqpWorkingSet* working_set);

/**
 * Computes the residuals of the given constraint values with respect
 * to the upper / lower bounds of the underlying problem.
 * The residuals are signed, i.e., positive if a constraint value exceeds its
 * upper bound, and negative if the lower bound exceeds the constraint value.
 **/
SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_signed_feasibility_residuals(SleqpProblem* problem,
                                   const SleqpVec* cons_val,
                                   SleqpVec* residuals,
                                   SleqpWorkingSet* working_set);

/**
 * Compute the maximum violation of the given constraint values, i.e.,
 * the \f$\| \cdot \|_{\infty}\f$-norm of the violation
 *
 * @param[in]  problem       Underlying problem
 * @param[in]  cons_val      Constraint values
 * @param[out] max_violation Total violation
 **/
SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_max_violation(SleqpProblem* problem,
                    SleqpVec* cons_val,
                    double* max_violation);

/**
 * Compute the total violation of the given constraint values, i.e.,
 * the \f$\| \cdot \|_{1}\f$-norm of the violation
 *
 * @param[in]  problem         Underlying problem
 * @param[in]  cons_val        Constraint values
 * @param[out] total_violation Total violation
 **/
SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_total_violation(SleqpProblem* problem,
                      SleqpVec* cons_val,
                      double* total_violation);

#endif /* SLEQP_FEAS_H */
