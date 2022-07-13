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
 * minimum norm by solving the problem
 *
 * \f[
 * \min \|x\|_2, \st A_W x = b_W
 * \f]
 *
 * @param[in]  aug_jac    The augmented Jacobian system
 * @param[in]  rhs        The right hand side \f$ b_W \f$
 * @param[out] sol        The solution \f$ x \f$
 *
 **/
SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_aug_jac_solve_min_norm(SleqpAugJac* aug_jac,
                             const SleqpVec* rhs,
                             SleqpVec* sol);

/**
 * Solves the following least-squares problem:
 *
 * \f[
 * \min_{x} \|y - x^{T} A_W\|_2
 * \f]
 *
 * @param[in]  aug_jac    The augmented Jacobian system
 * @param[in]  rhs        The target \f$ y \f$
 * @param[out] sol        The solution \f$ x \f$
 *
 **/
SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_aug_jac_solve_lsq(SleqpAugJac* aug_jac,
                        const SleqpVec* rhs,
                        SleqpVec* sol);

/**
 * Computes the projection of a given \f$ y \f$ onto the
 * null space of \f$ A_W by solving the problem.
 *
 * \f[
 * \min_{x} \|x - y\|_2, \st A_W x = 0
 * \f]
 *
 * @param[in]  aug_jac    The augmented Jacobian system
 * @param[in]  rhs        The target \f$ y \f$
 * @param[out] sol        The solution \f$ x \f$
 *
 **/
SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_aug_jac_project_nullspace(SleqpAugJac* aug_jac,
                                const SleqpVec* rhs,
                                SleqpVec* sol);

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
