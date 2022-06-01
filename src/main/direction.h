#ifndef SLEQP_DIRECTION_H
#define SLEQP_DIRECTION_H

#include "params.h"
#include "problem.h"
#include "pub_iterate.h"
#include "pub_types.h"
#include "types.h"

typedef struct SleqpDirection SleqpDirection;

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_direction_create(SleqpDirection** star,
                       SleqpProblem* problem,
                       SleqpParams* params);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_direction_reset(SleqpDirection* direction,
                      SleqpProblem* problem,
                      const SleqpIterate* iterate,
                      const SleqpVec* multipliers,
                      double* cache,
                      double zero_eps);

SLEQP_RETCODE
sleqp_direction_check(const SleqpDirection* direction,
                      SleqpProblem* problem,
                      const SleqpIterate* iterate,
                      const SleqpVec* multipliers,
                      double* cache,
                      double zero_eps,
                      bool* valid);

SleqpVec*
sleqp_direction_primal(const SleqpDirection* direction);

double*
sleqp_direction_obj_grad(SleqpDirection* direction);

SleqpVec*
sleqp_direction_cons_jac(const SleqpDirection* direction);

SleqpVec*
sleqp_direction_hess(const SleqpDirection* direction);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_direction_scale(SleqpDirection* direction, double factor);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_direction_set_zero(SleqpDirection* direction);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_direction_add_scaled(const SleqpDirection* first,
                           const SleqpDirection* second,
                           const double first_factor,
                           const double second_factor,
                           const double eps,
                           SleqpDirection* result);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_direction_copy(const SleqpDirection* source, SleqpDirection* target);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_direction_capture(SleqpDirection* direction);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_direction_release(SleqpDirection** star);

#endif /* SLEQP_DIRECTION_H */
