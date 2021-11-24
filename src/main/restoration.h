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
 * \min \: & f(x)                       \\
 * \text{s.t. } \: & l \leq c(x) \leq u \\
 * & l_x \leq x \leq u_x
 * \end{aligned}
 * \f]
 *
 * is given in terms of the residuum \f$ r : \mathbb{R}^{m + n} \to
 *\mathbb{R}^{m} \f$ \f[ r(x, s) := c(x) - s \f]
 *
 * as the problem
 *
 * \f[
 * \begin{aligned}
 * \min \: & 1/2 \| r(x, s) \|^2  \\
 * \text{s.t. } \: & l s \leq u   \\
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
                                    const SleqpSparseVec* primal,
                                    const SleqpSparseVec* cons_val,
                                    SleqpSparseVec* result);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_restoration_func_cons_val(SleqpFunc* restoration_func,
                                SleqpSparseVec** star);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_restoration_problem_restore(SleqpProblem* problem,
                                  const SleqpSparseVec* input,
                                  SleqpSparseVec* result);

#endif /* RESTORATION_H */
