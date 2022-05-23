#ifndef SLEQP_PARAMETRIC_H
#define SLEQP_PARAMETRIC_H

/**
 * @file solver.h
 * @brief A parametric trust region solver.
 *
 * This is an implementation of the algorithm laid out in
 * "An active-set algorithm for nonlinear programming using parametric linear
 *programming".
 *
 * We use the approximate parametric solve as suggested
 * by the authors in order to simplify the implementation
 *
 **/

#include "iterate.h"
#include "linesearch.h"
#include "merit.h"
#include "options.h"

#include "cauchy/cauchy.h"

typedef struct SleqpParametricSolver SleqpParametricSolver;

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_parametric_solver_create(SleqpParametricSolver** star,
                               SleqpProblem* problem,
                               SleqpParams* params,
                               SleqpOptions* options,
                               SleqpMerit* merit,
                               SleqpLineSearch* linesearch);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_parametric_solver_set_penalty(SleqpParametricSolver* solver,
                                    double penalty_parameter);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_parametric_solver_solve(SleqpParametricSolver* solver,
                              SleqpIterate* iterate,
                              SleqpCauchy* cauchy_data,
                              const SleqpVec* cauchy_step,
                              const SleqpVec* multipliers,
                              SleqpDirection* cauchy_direction,
                              double* trust_radius,
                              double* quadratic_merit_value);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_parametric_solver_capture(SleqpParametricSolver* solver);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_parametric_solver_release(SleqpParametricSolver** star);

#endif /* SLEQP_PARAMETRIC_H */
