#ifndef SLEQP_FACT_TYPES_H
#define SLEQP_FACT_TYPES_H

#include "sparse/mat.h"
#include "sparse/vec.h"

typedef struct SleqpFact SleqpFact;

typedef SLEQP_RETCODE (*SLEQP_FACT_SET_MATRIX)(void* fact_data,
                                               SleqpMat* matrix);

typedef SLEQP_RETCODE (*SLEQP_FACT_SOLVE)(void* fact_data, const SleqpVec* rhs);

typedef SLEQP_RETCODE (*SLEQP_FACT_SOLUTION)(void* fact_data,
                                             SleqpVec* sol,
                                             int begin,
                                             int end,
                                             double zero_eps);

typedef SLEQP_RETCODE (*SLEQP_FACT_CONDITION)(void* fact_data,
                                              double* condition);

typedef SLEQP_RETCODE (*SLEQP_FACT_FREE)(void** star);

typedef struct
{
  SLEQP_FACT_SET_MATRIX set_matrix;
  SLEQP_FACT_SOLVE solve;
  SLEQP_FACT_SOLUTION solution;
  SLEQP_FACT_CONDITION condition;
  SLEQP_FACT_FREE free;
} SleqpFactCallbacks;

#endif /* SLEQP_FACT_TYPES_H */
