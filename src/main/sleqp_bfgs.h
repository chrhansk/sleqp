#ifndef SLEQP_BFGS_H
#define SLEQP_BFGS_H

/**
 * @file sleqp_bfgs.h
 * @brief Defintion of BFGS method.
 **/

#include "sleqp_problem.h"

#include "sleqp_iterate.h"
#include "sparse/sleqp_sparse_vec.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpBFGSData SleqpBFGSData;

  SLEQP_RETCODE sleqp_bfgs_data_create(SleqpBFGSData** star,
                                       SleqpFunc* func,
                                       const SleqpParams* params,
                                       int num,
                                       bool damped);

  SLEQP_RETCODE sleqp_bfgs_data_push(SleqpBFGSData* data,
                                     SleqpIterate* old_iterate,
                                     SleqpIterate* new_iterate,
                                     SleqpSparseVec* multipliers);

  SLEQP_RETCODE sleqp_bfgs_data_hess_prod(SleqpBFGSData* data,
                                          const SleqpSparseVec* direction,
                                          SleqpSparseVec* product);

  SleqpFunc* sleqp_bfgs_get_func(SleqpBFGSData* data);

  SLEQP_RETCODE sleqp_bfgs_data_capture(SleqpBFGSData* data);

  SLEQP_RETCODE sleqp_bfgs_data_release(SleqpBFGSData** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_BFGS_H */
