#ifndef SLEQP_LSQR_H
#define SLEQP_LSQR_H

#include "lsqr_types.h"
#include "params.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpLSQRSolver SleqpLSQRSolver;

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_lsqr_solver_create(SleqpLSQRSolver** star,
                                         SleqpParams* params,
                                         int forward_dim,
                                         int adjoint_dim,
                                         SleqpLSQRCallbacks* callbacks,
                                         void* data);

  SLEQP_RETCODE sleqp_lsqr_solver_resize(SleqpLSQRSolver* solver,
                                         int forward_dim,
                                         int adjoint_dim);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_lsqr_solver_solve(SleqpLSQRSolver* solver,
                                        const SleqpSparseVec* rhs,
                                        double rel_tol,
                                        double trust_radius,
                                        SleqpSparseVec* sol);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_lsqr_solver_capture(SleqpLSQRSolver* solver);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_lsqr_solver_release(SleqpLSQRSolver** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_LSQR_H */
