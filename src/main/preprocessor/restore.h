#ifndef SLEQP_RESTORE_H
#define SLEQP_RESTORE_H

#include "problem.h"
#include "iterate.h"

#include "preprocessor/preprocessing_state.h"

typedef struct SleqpRestoration SleqpRestoration;

SLEQP_RETCODE sleqp_restoration_create(SleqpRestoration** star,
                                       SleqpPreprocessingState* preprocessing_state,
                                       SleqpProblem* transformed_problem,
                                       SleqpParams* params);

SLEQP_RETCODE sleqp_restoration_restore_iterate(SleqpRestoration* restoration,
                                                const SleqpIterate* source,
                                                SleqpIterate* target);

SLEQP_RETCODE sleqp_restoration_capture(SleqpRestoration* restoration);

SLEQP_RETCODE sleqp_restoration_release(SleqpRestoration** star);


#endif /* SLEQP_RESTORE_H */
