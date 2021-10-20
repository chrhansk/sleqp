#ifndef SLEQP_NEWTON_H
#define SLEQP_NEWTON_H

/**
 * @file newton.h
 * @brief Definition of functions used for the computation of Newton (aka EQP) steps.
 **/

#include "options.h"
#include "params.h"
#include "problem.h"
#include "iterate.h"
#include "timer.h"
#include "working_step.h"

#include "aug_jac/aug_jac.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpNewtonSolver SleqpNewtonSolver;

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_newton_solver_create(SleqpNewtonSolver** star,
                                           SleqpProblem* problem,
                                           SleqpWorkingStep* step,
                                           SleqpParams* params,
                                           SleqpOptions* options);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_newton_set_time_limit(SleqpNewtonSolver* data,
                                            double time_limit);

  SleqpTimer* sleqp_newton_get_timer(SleqpNewtonSolver* data);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_newton_set_iterate(SleqpNewtonSolver* data,
                                         SleqpIterate* iterate,
                                         SleqpAugJac* jacobian,
                                         double trust_radius,
                                         double penalty_parameter);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_newton_objective(SleqpNewtonSolver* data,
                                       const SleqpSparseVec* multipliers,
                                       const SleqpSparseVec* direction,
                                       double* objective);


  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_newton_add_violated_multipliers(SleqpNewtonSolver* data,
                                                      SleqpSparseVec* multipliers);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_newton_compute_step(SleqpNewtonSolver* data,
                                          SleqpSparseVec* multipliers,
                                          SleqpSparseVec* newton_step);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_newton_current_rayleigh(SleqpNewtonSolver* data,
                                              double* min_rayleigh,
                                              double* max_rayleigh);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_newton_solver_capture(SleqpNewtonSolver* data);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_newton_solver_release(SleqpNewtonSolver** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_NEWTON_H */
