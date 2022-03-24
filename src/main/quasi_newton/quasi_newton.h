#ifndef SLEQP_QUASI_NEWTON_H
#define SLEQP_QUASI_NEWTON_H

#include "quasi_newton_types.h"

#include "options.h"
#include "params.h"
#include "timer.h"

#include "bfgs.h"
#include "sr1.h"

typedef struct SleqpQuasiNewton SleqpQuasiNewton;

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_quasi_newton_create(SleqpQuasiNewton** star,
                          SleqpFunc* func,
                          SleqpQuasiNewtonCallbacks* callbacks,
                          void* quasi_newton_data);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_quasi_newton_push(SleqpQuasiNewton* quasi_newton,
                        const SleqpIterate* old_iterate,
                        const SleqpIterate* new_iterate,
                        const SleqpVec* multipliers);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_quasi_newton_reset(SleqpQuasiNewton* quasi_newton);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_quasi_newton_hess_prod(SleqpQuasiNewton* quasi_newton,
                             const SleqpVec* direction,
                             SleqpVec* product);

SleqpTimer*
sleqp_quasi_newton_update_timer(SleqpQuasiNewton* quasi_newton);

SleqpFunc*
sleqp_quasi_newton_get_func(SleqpQuasiNewton* quasi_newton);

SLEQP_NODISCARD SLEQP_RETCODE
sleqp_quasi_newton_capture(SleqpQuasiNewton* quasi_newton);

SLEQP_NODISCARD SLEQP_RETCODE
sleqp_quasi_newton_release(SleqpQuasiNewton** star);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_quasi_newton_create_default(SleqpQuasiNewton** star,
                                  SleqpFunc* func,
                                  SleqpParams* params,
                                  SleqpOptions* options);

#endif /* SLEQP_QUASI_NEWTON_H */
