#ifndef LSQR_H
#define LSQR_H

#include "problem.h"
#include "working_step.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpLSQRSolver SleqpLSQRSolver;

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_lsqr_solver_create(SleqpLSQRSolver** star,
                                         SleqpProblem* problem,
                                         SleqpWorkingStep* step,
                                         SleqpParams* params);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_lsqr_solver_compute_step(SleqpLSQRSolver* solver,
                                               SleqpSparseVec* lsqr_step);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_lsqr_solver_set_iterate(SleqpLSQRSolver* solver,
                                              SleqpIterate* iterate,
                                              SleqpAugJac* jacobian,
                                              double trust_radius,
                                              double penalty_parameter);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_lsqr_solver_capture(SleqpLSQRSolver* solver);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_lsqr_solver_release(SleqpLSQRSolver** star);

#ifdef __cplusplus
}
#endif

#endif /* LSQR_H */
