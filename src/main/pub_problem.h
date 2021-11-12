#ifndef SLEQP_PUB_PROBLEM_H
#define SLEQP_PUB_PROBLEM_H

#include "sleqp/export.h"

#include "sleqp/pub_func.h"
#include "sleqp/pub_params.h"
#include "sparse/pub_sparse_matrix.h"
#include "sparse/pub_sparse_vec.h"

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
 * where \f$ f : \mathbb{R}^{n} \to \mathbb{R} \f$, \f$ c : \mathbb{R}^{n} \to \mathbb{R}^{m} \f$
 * are functions, \f$ l, u \in \mathbb{R}^{m}, l \leq u \f$ are the constraint bounds, and
 * \f$ l_x, u_x \in \mathbb{R}^{n}, l_x \leq u_x \f$ are the variable bounds.
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
SLEQP_EXPORT SLEQP_NODISCARD
SLEQP_RETCODE sleqp_problem_create_simple(SleqpProblem** star,
                                          SleqpFunc* func,
                                          SleqpParams* params,
                                          const SleqpSparseVec* var_lb,
                                          const SleqpSparseVec* var_ub,
                                          const SleqpSparseVec* general_lb,
                                          const SleqpSparseVec* general_ub);

SLEQP_EXPORT SLEQP_NODISCARD
SLEQP_RETCODE sleqp_problem_create(SleqpProblem** star,
                                   SleqpFunc* func,
                                   SleqpParams* params,
                                   const SleqpSparseVec* var_lb,
                                   const SleqpSparseVec* var_ub,
                                   const SleqpSparseVec* genereal_lb,
                                   const SleqpSparseVec* genereal_ub,
                                   const SleqpSparseMatrix* linear_coeffs,
                                   const SleqpSparseVec* linear_lb,
                                   const SleqpSparseVec* linear_ub);

/**
 * Returns the total number of constraints (both general and linear) of the problem.
 **/
SLEQP_EXPORT int sleqp_problem_num_constraints(SleqpProblem* problem);

/**
 * Returns the total number of linear constraints of the problem.
 **/
SLEQP_EXPORT int sleqp_problem_num_linear_constraints(SleqpProblem* problem);

/**
 * Returns the total number of general constraints of the problem.
 **/
SLEQP_EXPORT int sleqp_problem_num_general_constraints(SleqpProblem* problem);

SLEQP_EXPORT SleqpFunc* sleqp_problem_func(SleqpProblem* problem);

SLEQP_EXPORT int sleqp_problem_num_variables(SleqpProblem* problem);

SLEQP_EXPORT SleqpSparseVec* sleqp_problem_var_lb(SleqpProblem* problem);

SLEQP_EXPORT SleqpSparseVec* sleqp_problem_var_ub(SleqpProblem* problem);

SLEQP_EXPORT SleqpSparseVec* sleqp_problem_general_lb(SleqpProblem* problem);

SLEQP_EXPORT SleqpSparseVec* sleqp_problem_general_ub(SleqpProblem* problem);

SLEQP_EXPORT SleqpSparseMatrix* sleqp_problem_linear_coeffs(SleqpProblem* problem);

SLEQP_EXPORT SleqpSparseVec* sleqp_problem_linear_lb(SleqpProblem* problem);

SLEQP_EXPORT SleqpSparseVec* sleqp_problem_linear_ub(SleqpProblem* problem);

SLEQP_EXPORT SleqpSparseVec* sleqp_problem_cons_lb(SleqpProblem* problem);

SLEQP_EXPORT SleqpSparseVec* sleqp_problem_cons_ub(SleqpProblem* problem);

SLEQP_EXPORT SLEQP_NODISCARD
SLEQP_RETCODE sleqp_problem_capture(SleqpProblem* problem);

SLEQP_EXPORT SLEQP_NODISCARD
SLEQP_RETCODE sleqp_problem_release(SleqpProblem** star);

/**
 * @}
 **/

#endif /* SLEQP_PUB_PROBLEM_H */
