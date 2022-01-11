#ifndef SLEQP_AUG_JAC_H
#define SLEQP_AUG_JAC_H

#include "aug_jac_types.h"
#include "timer.h"

typedef struct SleqpAugJac SleqpAugJac;

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_aug_jac_create(SleqpAugJac** star,
                     SleqpProblem* problem,
                     SleqpAugJacCallbacks* callbacks,
                     void* aug_jac_data);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_aug_jac_set_iterate(SleqpAugJac* aug_jac, SleqpIterate* iterate);

/**
 * Computes the solution of the system \f$ A_W x = b_W \f$ with
 * minimum norm.
 *
 * \f[
 * \min \|x\|_2, \text{s.t. } A_W x = b_W
 * \f]
 *
 * @param[in]  aug_jac    The augmented Jacobian system
 * @param[in]  rhs        The right hand side \f$ b_W \f$
 * @param[out] sol        The solution \f$ x \f$
 *
 **/
SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_aug_jac_min_norm_solution(SleqpAugJac* aug_jac,
                                const SleqpSparseVec* rhs,
                                SleqpSparseVec* sol);

/**
 * Computes the projection of the right hand side onto the
 * null space of the active constraints. If \f$ A_W \f$
 * is the jacobian of the working set, \f$ x_0 \f$ the
 * right hand side, then the output will contain a
 * vector \f$ x \f$ which solves
 *
 * \f[
 * \min \|x - x_0\|_2, \text{s.t. } A_W x = 0
 * \f]
 *
 * @param[in]  aug_jac     The augmented Jacobian system
 * @param[in]  rhs         The right hand side \f$ x_0 \f$
 * @param[out] primal_sol  The primal solution \f$ x \f$
 * @param[out] dual_sol    The dual solution
 *
 **/
SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_aug_jac_projection(SleqpAugJac* aug_jac,
                         const SleqpSparseVec* rhs,
                         SleqpSparseVec* primal_sol,
                         SleqpSparseVec* dual_sol);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_aug_jac_condition(SleqpAugJac* aug_jac, bool* exact, double* condition);

SleqpTimer*
sleqp_aug_jac_creation_timer(SleqpAugJac* aug_jac);

SleqpTimer*
sleqp_aug_jac_solution_timer(SleqpAugJac* aug_jac);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_aug_jac_capture(SleqpAugJac* jacobian);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_aug_jac_release(SleqpAugJac** star);

#endif /* SLEQP_AUG_JAC_H */
