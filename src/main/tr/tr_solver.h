#ifndef SLEQP_TR_SOLVER_H
#define SLEQP_TR_SOLVER_H

/**
 * @file tr_solver.h
 * @brief Definition of the EQP subproblem solver used to compute Newton (aka EQP) steps.
 **/

#include "options.h"
#include "params.h"
#include "problem.h"
#include "iterate.h"

#include "tr_types.h"

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
                                      SleqpAugJac* jacobian,
                                      SleqpSparseVec* multipliers,
                                      SleqpSparseVec* gradient,
                                      SleqpSparseVec* newton_step,
                                      double trust_radius,
                                      double* tr_dual);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_tr_solver_current_rayleigh(SleqpTRSolver* solver,
                                                 double* min_rayleigh,
                                                 double* max_rayleigh);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_tr_solver_capture(SleqpTRSolver* solver);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_tr_solver_release(SleqpTRSolver** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_TR_SOLVER_H */
