#ifndef SLEQP_EQP_TYPES_H
#define SLEQP_EQP_TYPES_H

#include "aug_jac/aug_jac.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef SLEQP_RETCODE (*SLEQP_EQP_SET_ITERATE)(SleqpIterate* iterate,
                                                 SleqpAugJac* jacobian,
                                                 double trust_radius,
                                                 double penalty_parameter,
                                                 void* eqp_data);

  typedef SLEQP_RETCODE (*SLEQP_EQP_SET_TIME_LIMIT)(double time_limit,
                                                    void* eqp_data);

  typedef SLEQP_RETCODE (*SLEQP_EQP_ADD_VIOLATED_MULTIPLIERS)(SleqpSparseVec* multipliers,
                                                              void* eqp_data);

  typedef SLEQP_RETCODE (*SLEQP_EQP_COMPUTE_STEP)(const SleqpSparseVec* multipliers,
                                                  SleqpSparseVec* newton_step,
                                                  void* eqp_data);

  typedef SLEQP_RETCODE (*SLEQP_EQP_CURRENT_RAYLEIGH)(double* min_rayleigh,
                                                      double* max_rayleigh,
                                                      void* eqp_data);

  typedef SLEQP_RETCODE (*SLEQP_EQP_FREE)(void* eqp_data);

  typedef struct {
    SLEQP_EQP_SET_ITERATE set_iterate;
    SLEQP_EQP_SET_TIME_LIMIT set_time_limit;
    SLEQP_EQP_ADD_VIOLATED_MULTIPLIERS add_violated_multipliers;
    SLEQP_EQP_COMPUTE_STEP compute_step;
    SLEQP_EQP_CURRENT_RAYLEIGH current_rayleigh;
    SLEQP_EQP_FREE free;

  } SleqpEQPCallbacks;

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_EQP_TYPES_H */
