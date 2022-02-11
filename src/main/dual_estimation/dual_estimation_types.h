#ifndef SLEQP_DUAL_ESTIMATION_TYPES_H
#define SLEQP_DUAL_ESTIMATION_TYPES_H

#include "iterate.h"
#include "working_set.h"

typedef SLEQP_RETCODE (*SLEQP_ESTIMATE_DUALS)(const SleqpIterate* iterate,
                                              SleqpSparseVec* cons_dual,
                                              SleqpSparseVec* vars_dual,
                                              void* estimation_data);

typedef SLEQP_RETCODE (*SLEQP_ESTIMATION_FREE)(void* estimation_data);

typedef struct
{
  SLEQP_ESTIMATE_DUALS estimate_duals;
  SLEQP_ESTIMATION_FREE estimation_free;
} SleqpDualEstimationCallbacks;

#endif /* SLEQP_DUAL_ESTIMATION_TYPES_H */
