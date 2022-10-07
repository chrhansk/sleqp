#ifndef SLEQP_TYPES_H
#define SLEQP_TYPES_H

#include "pub_types.h"

#include "enum.h"

typedef enum
{
  SLEQP_SOLVER_PHASE_OPTIMIZATION = 0,
  SLEQP_SOLVER_PHASE_RESTORATION,
  SLEQP_SOLVER_NUM_PHASES
} SLEQP_SOLVER_PHASE;

const SleqpEnum*
sleqp_enum_active_state();

const SleqpEnum*
sleqp_enum_status();

const SleqpEnum*
sleqp_enum_deriv_check();

const SleqpEnum*
sleqp_enum_hess_eval();

const SleqpEnum*
sleqp_enum_bfgs_sizing();

const SleqpEnum*
sleqp_enum_steptype();

const SleqpEnum*
sleqp_enum_dual_estimation();

const SleqpEnum*
sleqp_enum_tr_solver();

const SleqpEnum*
sleqp_enum_polishing_type();

const SleqpEnum*
sleqp_enum_parametric_cauchy();

const SleqpEnum*
sleqp_enum_initial_tr();

const SleqpEnum*
sleqp_enum_linesearch();

const SleqpEnum*
sleqp_enum_step_rule();

const SleqpEnum*
sleqp_enum_aug_jac_method();

#endif /* SLEQP_TYPES_H */
