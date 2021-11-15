#ifndef SLEQP_LSQR_H
#define SLEQP_LSQR_H

/**
 * @file lsqr.h
 * @brief Solve a linear least-squares problem.
 *
 * Solves \f$ \min_{x} \|Ax - b\|_2^2 \f$ where
 * \f$ A \f$ is given in terms of forward and
 * adjoint products using the LSQR iterative
 * algorithm. Optionally, a trust region can be
 * specified, in which case the iterations are
 * aborted as soon as the iterates leave the
 * trust region.
 **/

#include "lsqr_types.h"
#include "params.h"

typedef struct SleqpLSQRSolver SleqpLSQRSolver;

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_lsqr_solver_create(SleqpLSQRSolver** star,
                         SleqpParams* params,
                         int forward_dim,
                         int adjoint_dim,
                         SleqpLSQRCallbacks* callbacks,
                         void* data);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_lsqr_set_time_limit(SleqpLSQRSolver* solver, double time_limit);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_lsqr_solver_resize(SleqpLSQRSolver* solver,
                         int forward_dim,
                         int adjoint_dim);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_lsqr_solver_solve(SleqpLSQRSolver* solver,
                        const SleqpSparseVec* rhs,
                        double rel_tol,
                        double trust_radius,
                        SleqpSparseVec* sol);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_lsqr_solver_capture(SleqpLSQRSolver* solver);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_lsqr_solver_release(SleqpLSQRSolver** star);

#endif /* SLEQP_LSQR_H */
