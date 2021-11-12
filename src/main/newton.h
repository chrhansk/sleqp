#ifndef SLEQP_NEWTON_H
#define SLEQP_NEWTON_H

/**
 * @file newton.h
 * @brief Definition of functions used for the computation of Newton (aka EQP) steps.
 **/

#include "eqp.h"
#include "options.h"
#include "params.h"
#include "problem.h"
#include "iterate.h"
#include "timer.h"
#include "working_step.h"

#include "aug_jac/aug_jac.h"

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_newton_solver_create(SleqpEQPSolver** star,
                                         SleqpProblem* problem,
                                         SleqpParams* params,
                                         SleqpOptions* options,
                                         SleqpWorkingStep* step);

#endif /* SLEQP_NEWTON_H */
