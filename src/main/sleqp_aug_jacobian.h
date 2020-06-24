#ifndef SLEQP_AUG_JACOBIAN_H
#define SLEQP_AUG_JACOBIAN_H

/**
 * @file sleqp_aug_jacobian.h
 * @brief Definition of the augmented Jacobian matrix.
 **/

#include "sleqp_problem.h"
#include "sleqp_iterate.h"
#include "sleqp_params.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpAugJacobian SleqpAugJacobian;

  SLEQP_RETCODE sleqp_aug_jacobian_create(SleqpAugJacobian** star,
                                          SleqpProblem* problem,
                                          SleqpParams* params);

  SLEQP_RETCODE sleqp_aug_jacobian_free(SleqpAugJacobian** star);

  SLEQP_RETCODE sleqp_aug_jacobian_set_iterate(SleqpAugJacobian* jacobian,
                                               SleqpIterate* iterate);

  SleqpTimer* sleqp_aug_jacobian_get_factorization_timer(SleqpAugJacobian* jacobian);

  SLEQP_RETCODE sleqp_aug_jacobian_get_condition_estimate(SleqpAugJacobian* jacobian,
                                                          double* condition_estimate);

  /**
   * computes the solution of the system \f$ A_W x = b_W \f$ with
   * minimum norm.
   *
   * \f[
   * \min \|x\|_2, \text{s.t. } A_W x = b_W
   * \f]
   *
   * @param[in]  jacobian   The augmented Jacobian
   * @param[in]  rhs        The right hand side \f$ b_W \f$
   * @param[out] sol        The solution \f$ x \f$
   *
   **/
  SLEQP_RETCODE sleqp_aug_jacobian_min_norm_solution(SleqpAugJacobian* jacobian,
                                                     SleqpSparseVec* rhs,
                                                     SleqpSparseVec* sol);

  /**
   * Compute the projection of the right hand side onto the
   * null space of the active constraints. If \f$ A_W \f$
   * is the jacobian of the working set, \f$ x_0 \f$ the
   * right hand side, then the output will contain a
   * vector \f$ x \f$ which solves
   *
   * \f[
   * \min \|x - x_0\|_2, \text{s.t. } A_W x = 0
   * \f]
   *
   * @param[in]  jacobian    The augmented Jacobian
   * @param[in]  rhs         The right hand side \f$ x_0 \f$
   * @param[out] primal_sol  The primal solution \f$ x \f$
   * @param[out] dual_sol    The dual solution
   *
   **/
  SLEQP_RETCODE sleqp_aug_jacobian_projection(SleqpAugJacobian* jacobian,
                                              SleqpSparseVec* rhs,
                                              SleqpSparseVec* primal_sol,
                                              SleqpSparseVec* dual_sol);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_AUG_JACOBIAN_H */
