#ifndef SLEQP_TR_SOLVER_H
#define SLEQP_TR_SOLVER_H

/**
 * @file sleqp_tr_solver.h
 * @brief Definition of the EQP subproblem solver used to compute Newton (aka EQP) steps.
 **/

#include "sleqp_aug_jacobian.h"
#include "sleqp_options.h"
#include "sleqp_params.h"
#include "sleqp_problem.h"
#include "sleqp_iterate.h"

#include "sleqp_tr_types.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpTRSolver SleqpTRSolver;

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_tr_solver_create(SleqpTRSolver** star,
                                       SleqpTRCallbacks* callbacks,
                                       void* solver_data);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_tr_solver_set_time_limit(SleqpTRSolver* solver,
                                               double time_limit);

  SleqpTimer* sleqp_tr_solver_get_solve_timer(SleqpTRSolver* solver);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_tr_solver_solve(SleqpTRSolver* solver,
                                      SleqpAugJacobian* jacobian,
                                      SleqpSparseVec* multipliers,
                                      SleqpSparseVec* gradient,
                                      SleqpSparseVec* newton_step,
                                      double trust_radius);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_tr_solver_capture(SleqpTRSolver* solver);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_tr_solver_release(SleqpTRSolver** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_TR_SOLVER_H */
