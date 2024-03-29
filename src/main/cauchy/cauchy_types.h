#ifndef SLEQP_CAUCHY_TYPES
#define SLEQP_CAUCHY_TYPES

#include "types.h"

#include "iterate.h"

typedef enum
{
  SLEQP_CAUCHY_OBJTYPE_DEFAULT,
  SLEQP_CAUCHY_OBJTYPE_FEAS,
  SLEQP_CAUCHY_OBJTYPE_MIXED,
  SLEQP_CAUCHY_NUM_OBJTYPES
} SLEQP_CAUCHY_OBJTYPE;

typedef SLEQP_RETCODE (*SLEQP_CAUCHY_SET_ITERATE)(SleqpIterate* iterate,
                                                  double trust_radius,
                                                  void* cauchy_data);

typedef SLEQP_RETCODE (*SLEQP_CAUCHY_SET_TRUST_RADIUS)(double trust_radius,
                                                       void* cauchy_data);

typedef SLEQP_RETCODE (*SLEQP_CAUCHY_SOLVE)(
  SleqpVec* gradient,
  double penalty,
  SLEQP_CAUCHY_OBJTYPE objective_type,
  void* cauchy_data);

typedef SLEQP_RETCODE (*SLEQP_CAUCHY_OBJ_VAL)(double* objective_value,
                                              void* cauchy_data);

typedef SLEQP_RETCODE (*SLEQP_CAUCHY_WORKING_SET)(SleqpIterate* iterate,
                                                  void* cauchy_data);

typedef SLEQP_RETCODE (*SLEQP_CAUCHY_LP_STEP)(SleqpVec* direction,
                                              void* cauchy_data);

typedef SLEQP_RETCODE (*SLEQP_CAUCHY_LOCALLY_INFEASIBLE)(
  bool* locally_infeasible,
  void* cauchy_data);

typedef SLEQP_RETCODE (*SLEQP_CAUCHY_ESTIMATE_DUALS)(
  const SleqpWorkingSet* working_set,
  SleqpVec* cons_dual,
  SleqpVec* vars_dual,
  void* cauchy_data);

typedef SLEQP_RETCODE (*SLEQP_CAUCHY_SET_TIME_LIMIT)(double time_limit,
                                                     void* cauchy_data);

typedef SLEQP_RETCODE (*SLEQP_CAUCHY_VIOLATION)(double* violation,
                                                void* cauchy_data);

typedef SLEQP_RETCODE (*SLEQP_CAUCHY_BASIS_CONDITION)(bool* exact,
                                                      double* condition,
                                                      void* cauchy_data);

typedef SLEQP_RETCODE (*SLEQP_CAUCHY_PRINT_STATS)(double total_elapsed,
                                                  void* cauchy_data);

typedef SLEQP_RETCODE (*SLEQP_CAUCHY_FREE)(void* cauchy_data);

typedef struct
{
  SLEQP_CAUCHY_SET_ITERATE set_iterate;
  SLEQP_CAUCHY_SET_TRUST_RADIUS set_trust_radius;
  SLEQP_CAUCHY_SOLVE solve;
  SLEQP_CAUCHY_OBJ_VAL obj_val;
  SLEQP_CAUCHY_WORKING_SET working_set;
  SLEQP_CAUCHY_LP_STEP lp_step;
  SLEQP_CAUCHY_LOCALLY_INFEASIBLE locally_infeasible;
  SLEQP_CAUCHY_ESTIMATE_DUALS estimate_duals;
  SLEQP_CAUCHY_VIOLATION violation;
  SLEQP_CAUCHY_SET_TIME_LIMIT set_time_limit;
  SLEQP_CAUCHY_BASIS_CONDITION basis_condition;
  SLEQP_CAUCHY_PRINT_STATS print_stats;
  SLEQP_CAUCHY_FREE free;
} SleqpCauchyCallbacks;

#endif /* SLEQP_CAUCHY_TYPES */
