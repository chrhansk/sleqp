#ifndef SLEQP_ITERATE_H
#define SLEQP_ITERATE_H

/**
 * @file sleqp_iterate.h
 * @brief Definition of iterate.
 **/


#include "sleqp_export.h"
#include "sleqp_problem.h"
#include "sleqp_working_set.h"

#include "sparse/sleqp_sparse_matrix.h"
#include "sparse/sleqp_sparse_vec.h"

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
  SLEQP_EXPORT SleqpSparseVec* sleqp_iterate_get_primal(SleqpIterate* iterate);

  /**
   * The current function value
   **/
  SLEQP_EXPORT double sleqp_iterate_get_func_val(SleqpIterate* iterate);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_iterate_set_func_val(SleqpIterate* iterate,
                                           double value);

  /**
   * The current function gradient. Has dimension = num_variables.
   **/
  SLEQP_EXPORT SleqpSparseVec* sleqp_iterate_get_func_grad(SleqpIterate* iterate);

  /**
   * The current constraint values. Has dimension = num_constraints.
   **/
  SLEQP_EXPORT SleqpSparseVec* sleqp_iterate_get_cons_val(SleqpIterate* iterate);

  /**
   * The Jacobian of the constraitns at the current iterate.
   * Has num_constraints many rows, num_variables many columns.
   */
  SLEQP_EXPORT SleqpSparseMatrix* sleqp_iterate_get_cons_jac(SleqpIterate* iterate);

  /**
   * The current working set.
   **/
  SLEQP_EXPORT SleqpWorkingSet* sleqp_iterate_get_working_set(SleqpIterate* iterate);

  /**
   * The dual values of the constraints. Has dimension = num_constraints.
   */
  SLEQP_EXPORT SleqpSparseVec* sleqp_iterate_get_cons_dual(SleqpIterate* iterate);

  /**
   * The dual values of the variable bounds. Has dimension = num_variables.
   */
  SLEQP_EXPORT SleqpSparseVec* sleqp_iterate_get_vars_dual(SleqpIterate* iterate);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_iterate_slackness_residuum(SleqpIterate* iterate,
                                                 SleqpProblem* problem,
                                                 double* slackness_residuum);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_iterate_feasibility_residuum(SleqpIterate* iterate,
                                                   SleqpProblem* problem,
                                                   double feas_eps,
                                                   double* feasibility_residuum);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_iterate_stationarity_residuum(SleqpIterate* iterate,
                                                    SleqpProblem* problem,
                                                    double* cache,
                                                    double* stationarity_residuum);

  SLEQP_EXPORT bool sleqp_iterate_is_feasible(SleqpIterate* iterate,
                                              double feasibility_residuum,
                                              double feasibility_tolerance);

  SLEQP_EXPORT bool sleqp_iterate_is_optimal(SleqpIterate* iterate,
                                             SleqpParams* params,
                                             double feasibility_residuum,
                                             double slackness_residuum,
                                             double stationarity_residuum);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_iterate_get_violated_constraints(SleqpIterate* iterate,
                                                       SleqpProblem* problem,
                                                       int* violated_constraints,
                                                       int* num_violated_constraints,
                                                       double feas_eps);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_iterate_copy(SleqpIterate* source,
                                   SleqpIterate* target);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_iterate_capture(SleqpIterate* iterate);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_iterate_release(SleqpIterate** star);


#ifdef __cplusplus
}
#endif

#endif /* SLEQP_ITERATE_H */
