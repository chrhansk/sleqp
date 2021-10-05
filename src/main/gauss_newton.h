#ifndef SLEQP_GAUSS_NEWTON_H
#define SLEQP_GAUSS_NEWTON_H

#include "problem.h"
#include "working_step.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpGaussNewtonSolver SleqpGaussNewtonSolver;

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_gauss_newton_solver_create(SleqpGaussNewtonSolver** star,
                                                 SleqpProblem* problem,
                                                 SleqpWorkingStep* step,
                                                 SleqpParams* params);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_gauss_newton_solver_compute_step(SleqpGaussNewtonSolver* solver,
                                                       SleqpSparseVec* lsqr_step);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_gauss_newton_solver_set_iterate(SleqpGaussNewtonSolver* solver,
                                                      SleqpIterate* iterate,
                                                      SleqpAugJac* jacobian,
                                                      double trust_radius,
                                                      double penalty_parameter);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_gauss_newton_solver_capture(SleqpGaussNewtonSolver* solver);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_gauss_newton_solver_release(SleqpGaussNewtonSolver** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_GAUSS_NEWTON_H */
