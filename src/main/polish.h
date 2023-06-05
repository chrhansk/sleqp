#ifndef SLEQP_POLISH_H
#define SLEQP_POLISH_H

#include "iterate.h"

typedef struct SleqpPolishing SleqpPolishing;

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_polishing_create(SleqpPolishing** star,
                       SleqpProblem* problem,
                       SleqpSettings* settings);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_polishing_polish(SleqpPolishing* polishing,
                       SleqpIterate* iterate,
                       SLEQP_POLISHING_TYPE polishing_type);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_polishing_capture(SleqpPolishing* polishing);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_polishing_release(SleqpPolishing** star);

#endif /* SLEQP_POLISH_H */
