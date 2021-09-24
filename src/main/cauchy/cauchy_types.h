#ifndef SLEQP_CAUCHY_TYPES
#define SLEQP_CAUCHY_TYPES

#include "types.h"

#include "iterate.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef enum {
    SLEQP_CAUCHY_OBJECTIVE_TYPE_DEFAULT,
    SLEQP_CAUCHY_OBJECTIVE_TYPE_FEASIBILITY,
    SLEQP_CAUCHY_OBJECTIVE_TYPE_MIXED,
    SLEQP_NUM_CAUCHY_OBJECTIVES
  } SLEQP_CAUCHY_OBJECTIVE_TYPE;

  typedef SLEQP_RETCODE (*SLEQP_CAUCHY_SET_ITERATE)(SleqpIterate* iterate,
                                                    double trust_radius,
                                                    void* cauchy_data);

  typedef SLEQP_RETCODE (*SLEQP_CAUCHY_SET_TRUST_RADIUS)(double trust_radius,
                                                         void* cauchy_data);

  typedef SLEQP_RETCODE (*SLEQP_CAUCHY_SOLVE)(SleqpSparseVec* gradient,
                                              double penalty,
                                              SLEQP_CAUCHY_OBJECTIVE_TYPE objective_type,
                                              void* cauchy_data);

  typedef SLEQP_RETCODE (*SLEQP_CAUCHY_GET_OBJECTIVE_VALUE)(double* objective_value,
                                                            void* cauchy_data);

  typedef SLEQP_RETCODE (*SLEQP_CAUCHY_GET_WORKING_SET)(SleqpIterate* iterate,
                                                        void* cauchy_data);

  typedef SLEQP_RETCODE (*SLEQP_CAUCHY_GET_DIRECTION)(SleqpSparseVec* direction,
                                                      void* cauchy_data);

  typedef SLEQP_RETCODE (*SLEQP_CAUCHY_LOCALLY_INFEASIBLE)(bool* locally_infeasible,
                                                           void* cauchy_data);

  typedef SLEQP_RETCODE (*SLEQP_CAUCHY_GET_DUAL_ESTIMATION)(SleqpIterate* iterate,
                                                            void* cauchy_data);

  typedef SLEQP_RETCODE (*SLEQP_CAUCHY_GET_VIOLATION)(double* violation,
                                                      void* cauchy_data);

  typedef SLEQP_RETCODE (*SLEQP_CAUCHY_GET_BASIS_CONDITION)(bool* exact,
                                                            double* condition,
                                                            void* cauchy_data);

  typedef SLEQP_RETCODE (*SLEQP_CAUCHY_FREE)(void* cauchy_data);

  typedef struct {
    SLEQP_CAUCHY_SET_ITERATE set_iterate;
    SLEQP_CAUCHY_SET_TRUST_RADIUS set_trust_radius;
    SLEQP_CAUCHY_SOLVE solve;
    SLEQP_CAUCHY_GET_OBJECTIVE_VALUE get_objective_value;
    SLEQP_CAUCHY_GET_WORKING_SET get_working_set;
    SLEQP_CAUCHY_GET_DIRECTION get_direction;
    SLEQP_CAUCHY_LOCALLY_INFEASIBLE locally_infeasible;
    SLEQP_CAUCHY_GET_DUAL_ESTIMATION get_dual_estimation;
    SLEQP_CAUCHY_GET_VIOLATION get_violation;
    SLEQP_CAUCHY_GET_BASIS_CONDITION get_basis_condition;
    SLEQP_CAUCHY_FREE free;
  } SleqpCauchyCallbacks;

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_CAUCHY_TYPES */
