#ifndef SLEQP_PROBLEM_H
#define SLEQP_PROBLEM_H

/**
 * @file sleqp_problem.h
 * @brief Definition of the programming problem.
 **/

#include "sleqp_func.h"
#include "sleqp_params.h"
#include "sleqp_types.h"

#include "sparse/sleqp_sparse_matrix.h"
#include "sparse/sleqp_sparse_vec.h"

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * @defgroup problem Problem definition
   *
   * A NLP is given as
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

  SleqpFunc* sleqp_problem_func(SleqpProblem* problem);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_problem_set_value(SleqpProblem* problem,
                                        SleqpSparseVec* x,
                                        SLEQP_VALUE_REASON reason,
                                        int* func_grad_nnz,
                                        int* cons_val_nnz,
                                        int* cons_jac_nnz);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_problem_eval(SleqpProblem* problem,
                                   const SleqpSparseVec* cons_indices,
                                   double* func_val,
                                   SleqpSparseVec* func_grad,
                                   SleqpSparseVec* cons_val,
                                   SleqpSparseMatrix* cons_jac);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_problem_val(SleqpProblem* problem,
                                  double* func_val);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_problem_grad(SleqpProblem* problem,
                                   SleqpSparseVec* func_grad);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_problem_cons_val(SleqpProblem* problem,
                                       const SleqpSparseVec* cons_indices,
                                       SleqpSparseVec* cons_val);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_problem_cons_jac(SleqpProblem* problem,
                                       const SleqpSparseVec* cons_indices,
                                       SleqpSparseMatrix* cons_jac);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_problem_hess_prod(SleqpProblem* problem,
                                        const double* func_dual,
                                        const SleqpSparseVec* direction,
                                        const SleqpSparseVec* cons_duals,
                                        SleqpSparseVec* product);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_problem_hess_bilinear(SleqpProblem* problem,
                                            const double* func_dual,
                                            const SleqpSparseVec* direction,
                                            const SleqpSparseVec* cons_duals,
                                            double* bilinear_prod);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_problem_capture(SleqpProblem* problem);
  
  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_problem_release(SleqpProblem** star);

  /**
   * @}
   **/

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_PROBLEM_H */
