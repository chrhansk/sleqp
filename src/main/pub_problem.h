#ifndef SLEQP_PUB_PROBLEM_H
#define SLEQP_PUB_PROBLEM_H

#include "sleqp/export.h"

#include "sleqp/pub_func.h"
#include "sleqp/pub_params.h"
#include "sparse/pub_sparse_matrix.h"
#include "sparse/pub_vec.h"

/**
 * @file pub_problem.h
 * @brief Definition of the programming problem.
 **/

/**
 * @defgroup problem Problem definition
 *
 * An NLP is given as
 *
 * \f[
 * \begin{aligned}
 * \min \: & f(x)                       \\
 * \text{s.t. } \: & l \leq c(x) \leq u \\
 * & l_x \leq x \leq u_x
 * \end{aligned}
 * \f]
 *
 * where \f$ f : \mathbb{R}^{n} \to \mathbb{R} \f$, \f$ c : \mathbb{R}^{n} \to
 *\mathbb{R}^{m} \f$ are functions, \f$ l, u \in \mathbb{R}^{m}, l \leq u \f$
 *are the constraint bounds, and \f$ l_x, u_x \in \mathbb{R}^{n}, l_x \leq u_x
 *\f$ are the variable bounds.
 *
 * @see Functions
 *
 * @{
 **/

typedef struct SleqpProblem SleqpProblem;

/**
 * Creates a new problem, which will be saved
 * in the newly allocated pointer.
 *
 * The function pointer is kept as a reference,
 * The upper / lower bound vectors are copied.
 *
 **/
SLEQP_EXPORT SLEQP_NODISCARD SLEQP_RETCODE
sleqp_problem_create_simple(SleqpProblem** star,
                            SleqpFunc* func,
                            SleqpParams* params,
                            const SleqpVec* var_lb,
                            const SleqpVec* var_ub,
                            const SleqpVec* general_lb,
                            const SleqpVec* general_ub);

SLEQP_EXPORT SLEQP_NODISCARD SLEQP_RETCODE
sleqp_problem_create(SleqpProblem** star,
                     SleqpFunc* func,
                     SleqpParams* params,
                     const SleqpVec* var_lb,
                     const SleqpVec* var_ub,
                     const SleqpVec* genereal_lb,
                     const SleqpVec* genereal_ub,
                     const SleqpSparseMatrix* linear_coeffs,
                     const SleqpVec* linear_lb,
                     const SleqpVec* linear_ub);

/**
 * Returns the total number of constraints (both general and linear) of the
 *problem.
 **/
SLEQP_EXPORT int
sleqp_problem_num_cons(SleqpProblem* problem);

/**
 * Returns the total number of linear constraints of the problem.
 **/
SLEQP_EXPORT int
sleqp_problem_num_lin_cons(SleqpProblem* problem);

/**
 * Returns the total number of general constraints of the problem.
 **/
SLEQP_EXPORT int
sleqp_problem_num_gen_cons(SleqpProblem* problem);

SLEQP_EXPORT SleqpFunc*
sleqp_problem_func(SleqpProblem* problem);

SLEQP_EXPORT int
sleqp_problem_num_vars(SleqpProblem* problem);

SLEQP_EXPORT SleqpVec*
sleqp_problem_vars_lb(SleqpProblem* problem);

SLEQP_EXPORT SleqpVec*
sleqp_problem_vars_ub(SleqpProblem* problem);

SLEQP_EXPORT SleqpVec*
sleqp_problem_general_lb(SleqpProblem* problem);

SLEQP_EXPORT SleqpVec*
sleqp_problem_general_ub(SleqpProblem* problem);

SLEQP_EXPORT SleqpSparseMatrix*
sleqp_problem_linear_coeffs(SleqpProblem* problem);

SLEQP_EXPORT SleqpVec*
sleqp_problem_linear_lb(SleqpProblem* problem);

SLEQP_EXPORT SleqpVec*
sleqp_problem_linear_ub(SleqpProblem* problem);

SLEQP_EXPORT SleqpVec*
sleqp_problem_cons_lb(SleqpProblem* problem);

SLEQP_EXPORT SleqpVec*
sleqp_problem_cons_ub(SleqpProblem* problem);

SLEQP_EXPORT SLEQP_NODISCARD SLEQP_RETCODE
sleqp_problem_capture(SleqpProblem* problem);

SLEQP_EXPORT SLEQP_NODISCARD SLEQP_RETCODE
sleqp_problem_release(SleqpProblem** star);

/**
 * @}
 **/

#endif /* SLEQP_PUB_PROBLEM_H */
