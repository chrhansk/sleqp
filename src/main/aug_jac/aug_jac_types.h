#ifndef SLEQP_AUG_JAC_TYPES_H
#define SLEQP_AUG_JAC_TYPES_H

#include "iterate.h"

typedef SLEQP_RETCODE (*SLEQP_AUG_JAC_SET_ITERATE)(SleqpIterate* iterate,
                                                   void* aug_jac);

typedef SLEQP_RETCODE (*SLEQP_AUG_JAC_SOLVE_MIN_NORM)(const SleqpVec* rhs,
                                                      SleqpVec* sol,
                                                      void* aug_jac);

typedef SLEQP_RETCODE (*SLEQP_AUG_JAC_SOLVE_LSQ)(const SleqpVec* rhs,
                                                 SleqpVec* sol,
                                                 void* aug_jac);

typedef SLEQP_RETCODE (*SLEQP_AUG_JAC_PROJECT_NULLSPACE)(const SleqpVec* rhs,
                                                         SleqpVec* sol,
                                                         void* aug_jac);

typedef SLEQP_RETCODE (*SLEQP_AUG_JAC_CONDITION)(bool* exact,
                                                 double* condition,
                                                 void* aug_jac);

typedef SLEQP_RETCODE (*SLEQP_AUG_JAC_FREE)(void* aug_jac);

typedef struct
{
  SLEQP_AUG_JAC_SET_ITERATE set_iterate;
  SLEQP_AUG_JAC_SOLVE_MIN_NORM solve_min_norm;
  SLEQP_AUG_JAC_SOLVE_LSQ solve_lsq;
  SLEQP_AUG_JAC_PROJECT_NULLSPACE project_nullspace;
  SLEQP_AUG_JAC_CONDITION condition;
  SLEQP_AUG_JAC_FREE free;
} SleqpAugJacCallbacks;

#endif /* SLEQP_AUG_JAC_TYPES_H */
