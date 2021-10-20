#ifndef SLEQP_GAUSS_NEWTON_H
#define SLEQP_GAUSS_NEWTON_H

#include "eqp.h"
#include "problem.h"
#include "working_step.h"

#ifdef __cplusplus
extern "C" {
#endif

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_gauss_newton_solver_create(SleqpEQPSolver** star,
                                                 SleqpProblem* problem,
                                                 SleqpParams* params,
                                                 SleqpWorkingStep* step);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_GAUSS_NEWTON_H */
