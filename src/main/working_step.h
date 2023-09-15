#ifndef SLEQP_WORKING_STEP_H
#define SLEQP_WORKING_STEP_H

/**
 * @file working_step.h
 * @brief Definition of functions used for the computation of the initial step
 *towards the working set.
 **/

#include "direction.h"
#include "iterate.h"
#include "problem.h"
#include "settings.h"

#include "aug_jac/aug_jac.h"

typedef struct SleqpWorkingStep SleqpWorkingStep;

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_working_step_create(SleqpWorkingStep** star,
                          SleqpProblem* problem,
                          SleqpSettings* settings);

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_working_step_set_iterate(SleqpWorkingStep* step,
                               SleqpIterate* iterate,
                               SleqpAugJac* jacobian,
                               double trust_radius);

SleqpVec*
sleqp_working_step_get_step(SleqpWorkingStep* step);

SleqpDirection*
sleqp_working_step_direction(SleqpWorkingStep* step);

double
sleqp_working_step_reduced_trust_radius(SleqpWorkingStep* step);

bool
sleqp_working_step_in_working_set(SleqpWorkingStep* step);

double
sleqp_working_step_newton_obj_offset(SleqpWorkingStep* step,
                                     double penalty_parameter);

SleqpVec*
sleqp_working_step_violated_cons_multipliers(SleqpWorkingStep* step);

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_working_step_set_multipliers(SleqpWorkingStep* step,
                                   const SleqpVec* multipliers);

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_working_step_capture(SleqpWorkingStep* step);

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_working_step_release(SleqpWorkingStep** star);

#endif /* SLEQP_WORKING_STEP_H */
