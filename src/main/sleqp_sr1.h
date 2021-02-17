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

  typedef struct SleqpSR1Data SleqpSR1Data;

  SLEQP_RETCODE sleqp_sr1_data_create(SleqpSR1Data** star,
                                      SleqpFunc* func,
                                      SleqpParams* params,
                                      SleqpOptions* options);

  SLEQP_RETCODE sleqp_sr1_data_push(SleqpSR1Data* data,
                                    SleqpIterate* old_iterate,
                                    SleqpIterate* new_iterate,
                                    SleqpSparseVec* multipliers);

  SLEQP_RETCODE sleqp_sr1_data_hess_prod(SleqpSR1Data* data,
                                         const SleqpSparseVec* direction,
                                         SleqpSparseVec* product);

  SleqpTimer* sleqp_sr1_update_timer(SleqpSR1Data* data);

  SleqpFunc* sleqp_sr1_get_func(SleqpSR1Data* data);

  SLEQP_RETCODE sleqp_sr1_data_capture(SleqpSR1Data* data);

  SLEQP_RETCODE sleqp_sr1_data_release(SleqpSR1Data** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_SR1_H */
