#ifndef SLEQP_TR_TYPES_H
#define SLEQP_TR_TYPES_H

#include "problem.h"

#include "aug_jac/aug_jac.h"
#include "sparse/sparse_vec.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef SLEQP_RETCODE (*SLEQP_TR_SOLVER_SOLVE)(SleqpAugJac* jacobian,
                                                 SleqpSparseVec* multipliers,
                                                 SleqpSparseVec* gradient,
                                                 SleqpSparseVec* newton_step,
                                                 double trust_radius,
                                                 double* tr_dual,
                                                 double time_limit,
                                                 void* solver_data);

  typedef SLEQP_RETCODE (*SLEQP_TR_SOLVER_RAYLEIGH)(double* min_rayleigh,
                                                    double* max_rayleigh,
                                                    void* solver_data);

  typedef SLEQP_RETCODE (*SLEQP_TR_SOLVER_FREE)(void** solver_data);

  typedef struct {
    SLEQP_TR_SOLVER_SOLVE solve;
    SLEQP_TR_SOLVER_RAYLEIGH rayleigh;
    SLEQP_TR_SOLVER_FREE free;
  } SleqpTRCallbacks;

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_TR_TYPES_H */
