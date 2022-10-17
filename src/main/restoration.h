#ifndef RESTORATION_H
#define RESTORATION_H

/**
 * @file restoration.h
 * @brief Defintion of the restoration problem.
 *
 * The restoration problem of the original problem
 *
 * \f[
 * \begin{aligned}
 * \min \: & f(x)              \\
 * \st \: & l \leq c(x) \leq u \\
 * & l_x \leq x \leq u_x
 * \end{aligned}
 * \f]
 *
 * is given in terms of the residuum \f$ r : \R^{m + n} \to
 *\R^{m} \f$ \f[ r(x, s) := c(x) - s \f]
 *
 * as the problem
 *
 * \f[
 * \begin{aligned}
 * \min \: & 1/2 \| r(x, s) \|^2 \\
 * \st \: & l s \leq u           \\
 * & l_x \leq x \leq u_x
 * \end{aligned}
 * \f]
 *
 **/

#include "iterate.h"
#include "problem.h"

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_restoration_problem_create(SleqpProblem** star,
                                 SleqpParams* params,
                                 SleqpProblem* problem);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_restoration_problem_transform(SleqpProblem* problem,
                                    const SleqpVec* primal,
                                    const SleqpVec* cons_val,
                                    SleqpVec* result);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_restoration_func_cons_val(SleqpFunc* restoration_func, SleqpVec** star);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_restoration_func_cons_jac(SleqpFunc* restoration_func,
                                SleqpSparseMatrix** star);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_restoration_func_init(SleqpFunc* restoration_func,
                            SleqpVec* restoration_primal,
                            SleqpVec* orig_cons_val,
                            SleqpSparseMatrix* orig_cons_jac);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_restoration_problem_restore(SleqpProblem* problem,
                                  const SleqpVec* input,
                                  SleqpVec* result);

#endif /* RESTORATION_H */
