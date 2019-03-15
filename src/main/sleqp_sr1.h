#ifndef SLEQP_SR1_H
#define SLEQP_SR1_H

#include "sleqp_func.h"

#include "sleqp_iterate.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpSR1Data SleqpSR1Data;

  SLEQP_RETCODE sleqp_sr1_data_create(SleqpSR1Data** star,
                                      SleqpFunc* func,
                                      SleqpParams* params,
                                      int num);

  SLEQP_RETCODE sleqp_sr1_data_push(SleqpSR1Data* data,
                                    SleqpIterate* old_iterate,
                                    SleqpIterate* new_iterate);

  SLEQP_RETCODE sleqp_sr1_data_hess_prod(SleqpSR1Data* data,
                                         SleqpSparseVec* direction,
                                         SleqpSparseVec* product);

  SleqpFunc* sleqp_sr1_get_func(SleqpSR1Data* data);

  SLEQP_RETCODE sleqp_sr1_data_free(SleqpSR1Data** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_SR1_H */
