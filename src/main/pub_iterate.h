#ifndef SLEQP_PUB_ITERATE_H
#define SLEQP_PUB_ITERATE_H

/**
 * @file pub_iterate.h
 * @brief Definition of iterate.
 **/


#include "sleqp/export.h"
#include "sleqp/pub_problem.h"
#include "sleqp/pub_working_set.h"

#include "sparse/pub_sparse_matrix.h"
#include "sparse/pub_sparse_vec.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpIterate SleqpIterate;

  /**
   * Create a new iterate
   *
   * @param[in,out] star    A pointer to the newly created iterate
   * @param[in] problem     The underlying problem
   * @param[in] x           The point of the iterate
   **/
  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_iterate_create(SleqpIterate** star,
                                     SleqpProblem* problem,
                                     SleqpSparseVec* x);

  /**
   * The current point. Has dimension = num_variables.
   **/
  SLEQP_EXPORT SleqpSparseVec* sleqp_iterate_get_primal(const SleqpIterate* iterate);

  /**
   * The current function value
   **/
  SLEQP_EXPORT double sleqp_iterate_get_func_val(const SleqpIterate* iterate);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_iterate_set_func_val(SleqpIterate* iterate,
                                           double value);

  /**
   * The current function gradient. Has dimension = num_variables.
   **/
  SLEQP_EXPORT SleqpSparseVec* sleqp_iterate_get_func_grad(const SleqpIterate* iterate);

  /**
   * The current constraint values. Has dimension = num_constraints.
   **/
  SLEQP_EXPORT SleqpSparseVec* sleqp_iterate_get_cons_val(const SleqpIterate* iterate);

  /**
   * The Jacobian of the constraitns at the current iterate.
   * Has num_constraints many rows, num_variables many columns.
   */
  SLEQP_EXPORT SleqpSparseMatrix* sleqp_iterate_get_cons_jac(const SleqpIterate* iterate);

  /**
   * The current working set.
   **/
  SLEQP_EXPORT SleqpWorkingSet* sleqp_iterate_get_working_set(const SleqpIterate* iterate);

  /**
   * The dual values of the constraints. Has dimension = num_constraints.
   */
  SLEQP_EXPORT SleqpSparseVec* sleqp_iterate_get_cons_dual(const SleqpIterate* iterate);

  /**
   * The dual values of the variable bounds. Has dimension = num_variables.
   */
  SLEQP_EXPORT SleqpSparseVec* sleqp_iterate_get_vars_dual(const SleqpIterate* iterate);


  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_iterate_capture(SleqpIterate* iterate);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_iterate_release(SleqpIterate** star);


#ifdef __cplusplus
}
#endif

#endif /* SLEQP_PUB_ITERATE_H */
