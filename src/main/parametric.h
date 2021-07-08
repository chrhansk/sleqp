#ifndef SLEQP_PARAMETRIC_H
#define SLEQP_PARAMETRIC_H

/**
 * @file solver.h
 * @brief A parametric trust region solver.
 *
 * This is an implementation of the algorithm laid out in
 * "An active-set algorithm for nonlinear programming using parametric linear programming".
 *
 * We use the approximate parametric solve as suggested
 * by the authors in order to simplify the implementation
 *
 **/

#include "cauchy.h"
#include "iterate.h"
#include "linesearch.h"
#include "merit.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpParametricSolver SleqpParametricSolver;

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_parametric_solver_create(SleqpParametricSolver** star,
                                               SleqpProblem* problem,
                                               SleqpParams* params,
                                               SleqpOptions* options,
                                               SleqpMeritData* merit_data,
                                               SleqpLineSearchData* linesearch);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_parametric_solver_set_penalty(SleqpParametricSolver* solver,
                                                    double penalty_parameter);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_parametric_solver_solve(SleqpParametricSolver* solver,
                                              SleqpIterate* iterate,
                                              SleqpCauchy* cauchy_data,
                                              SleqpSparseVec* cauchy_direction,
                                              SleqpSparseVec* cauchy_hessian_direction,
                                              const SleqpSparseVec* multipliers,
                                              double* trust_radius,
                                              double* quadratic_merit_value);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_parametric_solver_capture(SleqpParametricSolver* solver);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_parametric_solver_release(SleqpParametricSolver** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_PARAMETRIC_H */
