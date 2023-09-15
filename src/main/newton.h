#ifndef SLEQP_NEWTON_H
#define SLEQP_NEWTON_H

/**
 * @file newton.h
 * @brief Definition of functions used for the computation of Newton (aka EQP)
 *steps.
 **/

#include "eqp.h"
#include "iterate.h"
#include "problem.h"
#include "settings.h"
#include "timer.h"
#include "working_step.h"

#include "aug_jac/aug_jac.h"

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_newton_solver_create(SleqpEQPSolver** star,
                           SleqpProblem* problem,
                           SleqpSettings* settings,
                           SleqpWorkingStep* step);

#endif /* SLEQP_NEWTON_H */
