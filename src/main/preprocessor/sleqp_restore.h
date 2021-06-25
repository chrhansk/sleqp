#ifndef SLEQP_RESTORE_H
#define SLEQP_RESTORE_H

#include "sleqp_problem.h"
#include "sleqp_iterate.h"

#include "preprocessor/sleqp_preprocessing_state.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpRestoration SleqpRestoration;

  SLEQP_RETCODE sleqp_restoration_create(SleqpRestoration** star,
                                         SleqpPreprocessingState* preprocessing_state);

  SLEQP_RETCODE sleqp_restoration_restore_iterate(SleqpRestoration* restoration,
                                                  const SleqpIterate* source,
                                                  SleqpIterate* target);

  SLEQP_RETCODE sleqp_restoration_capture(SleqpRestoration* restoration);

  SLEQP_RETCODE sleqp_restoration_release(SleqpRestoration** star);


#ifdef __cplusplus
}
#endif

#endif /* SLEQP_RESTORE_H */
