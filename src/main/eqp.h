#ifndef SLEQP_EQP_H
#define SLEQP_EQP_H

#include "direction.h"
#include "eqp_types.h"
#include "timer.h"

typedef struct SleqpEQPSolver SleqpEQPSolver;

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_eqp_solver_create(SleqpEQPSolver** star,
                        SleqpEQPCallbacks* callbacks,
                        void* eqp_data);

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_eqp_solver_set_iterate(SleqpEQPSolver* solver,
                             SleqpIterate* iterate,
                             SleqpAugJac* jacobian,
                             double trust_radius,
                             double penalty_parameter);

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_eqp_solver_set_time_limit(SleqpEQPSolver* solver, double time_limit);

SleqpTimer*
sleqp_eqp_solver_get_timer(SleqpEQPSolver* solver);

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_eqp_solver_add_violated_multipliers(SleqpEQPSolver* solver,
                                          SleqpVec* multipliers);

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_eqp_solver_compute_direction(SleqpEQPSolver* solver,
                                   const SleqpVec* multipliers,
                                   SleqpDirection* newton_direction);

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_eqp_solver_current_rayleigh(SleqpEQPSolver* solver,
                                  double* min_rayleigh,
                                  double* max_rayleigh);

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_eqp_solver_capture(SleqpEQPSolver* solver);

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_eqp_solver_release(SleqpEQPSolver** star);

#endif /* SLEQP_EQP_H */
