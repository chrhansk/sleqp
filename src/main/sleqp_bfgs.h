#ifndef SLEQP_BFGS_H
#define SLEQP_BFGS_H

#include "sleqp_problem.h"

#include "sleqp_iterate.h"
#include "sparse/sleqp_sparse_vec.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpBFGSData SleqpBFGSData;

  SLEQP_RETCODE sleqp_bfgs_data_create(SleqpBFGSData** star,
                                       SleqpFunc* func,
                                       SleqpParams* params,
                                       int num,
                                       bool damped);

  SLEQP_RETCODE sleqp_bfgs_data_push(SleqpBFGSData* data,
                                     SleqpIterate* old_iterate,
                                     SleqpIterate* new_iterate);

  SLEQP_RETCODE sleqp_bfgs_data_hess_prod(SleqpBFGSData* data,
                                          SleqpSparseVec* direction,
                                          SleqpSparseVec* product);

  SleqpFunc* sleqp_bfgs_get_func(SleqpBFGSData* data);

  SLEQP_RETCODE sleqp_bfgs_data_free(SleqpBFGSData** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_BFGS_H */
