#ifndef SLEQP_FACTORIZATION_TYPES_H
#define SLEQP_FACTORIZATION_TYPES_H

#include "sparse/sparse_matrix.h"
#include "sparse/vec.h"

typedef struct SleqpFact SleqpFact;

typedef SLEQP_RETCODE (*SLEQP_FACT_SET_MATRIX)(void* factorization_data,
                                               SleqpSparseMatrix* matrix);

typedef SLEQP_RETCODE (*SLEQP_FACT_SOLVE)(void* factorization_data,
                                          const SleqpVec* rhs);

typedef SLEQP_RETCODE (*SLEQP_FACT_SOLUTION)(void* factorization_data,
                                             SleqpVec* sol,
                                             int begin,
                                             int end,
                                             double zero_eps);

typedef SLEQP_RETCODE (*SLEQP_FACT_CONDITION_ESTIMATE)(
  void* factorization_data,
  double* condition_estimate);

typedef SLEQP_RETCODE (*SLEQP_FACT_FREE)(void** star);

typedef struct
{
  SLEQP_FACT_SET_MATRIX set_matrix;
  SLEQP_FACT_SOLVE solve;
  SLEQP_FACT_SOLUTION solution;
  SLEQP_FACT_CONDITION_ESTIMATE condition_estimate;
  SLEQP_FACT_FREE free;
} SleqpFactorizationCallbacks;

#endif /* SLEQP_FACT                                                           \
_TYPES_H */
