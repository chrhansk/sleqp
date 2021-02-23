#ifndef SLEQP_TR_TYPES_H
#define SLEQP_TR_TYPES_H

#include "sleqp_aug_jacobian.h"
#include "sleqp_problem.h"

#include "sparse/sleqp_sparse_vec.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef SLEQP_RETCODE (*SLEQP_TR_SOLVER_SOLVE)(SleqpAugJacobian* jacobian,
                                                 SleqpSparseVec* multipliers,
                                                 SleqpSparseVec* gradient,
                                                 SleqpSparseVec* newton_step,
                                                 double trust_radius,
                                                 double time_limit,
                                                 void* solver_data);

  typedef SLEQP_RETCODE (*SLEQP_TR_SOLVER_FREE)(void** solver_data);

  typedef struct {
    SLEQP_TR_SOLVER_SOLVE solve;
    SLEQP_TR_SOLVER_FREE free;
  } SleqpTRCallbacks;

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_TR_TYPES_H */
