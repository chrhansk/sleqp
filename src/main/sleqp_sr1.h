#ifndef SLEQP_SR1_H
#define SLEQP_SR1_H

/**
 * @file sleqp_sr1.h
 * @brief Defintion of SR1 method.
 **/

#include "sleqp_func.h"
#include "sleqp_iterate.h"
#include "sleqp_options.h"
#include "sleqp_params.h"
#include "sleqp_timer.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpSR1 SleqpSR1;

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_sr1_create(SleqpSR1** star,
                                 SleqpFunc* func,
                                 SleqpParams* params,
                                 SleqpOptions* options);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_sr1_push(SleqpSR1* sr1,
                               SleqpIterate* old_iterate,
                               SleqpIterate* new_iterate,
                               SleqpSparseVec* multipliers);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_sr1_hess_prod(SleqpSR1* sr1,
                                    const SleqpSparseVec* direction,
                                    SleqpSparseVec* product);

  SleqpTimer* sleqp_sr1_update_timer(SleqpSR1* sr1);

  SleqpFunc* sleqp_sr1_get_func(SleqpSR1* sr1);

  SLEQP_NODISCARD SLEQP_RETCODE sleqp_sr1_capture(SleqpSR1* sr1);

  SLEQP_NODISCARD SLEQP_RETCODE sleqp_sr1_release(SleqpSR1** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_SR1_H */
