#ifndef SLEQP_FACTORIZATION_TYPES_H
#define SLEQP_FACTORIZATION_TYPES_H

#include "sparse/sparse_matrix.h"
#include "sparse/sparse_vec.h"

typedef struct SleqpFactorization SleqpFactorization;

typedef SLEQP_RETCODE (*SLEQP_FACTORIZATION_SET_MATRIX)(
  void* factorization_data,
  SleqpSparseMatrix* matrix);

typedef SLEQP_RETCODE (*SLEQP_FACTORIZATION_SOLVE)(void* factorization_data,
                                                   const SleqpSparseVec* rhs);

typedef SLEQP_RETCODE (*SLEQP_FACTORIZATION_SOLUTION)(void* factorization_data,
                                                      SleqpSparseVec* sol,
                                                      int begin,
                                                      int end,
                                                      double zero_eps);

typedef SLEQP_RETCODE (*SLEQP_FACTORIZATION_CONDITION_ESTIMATE)(
  void* factorization_data,
  double* condition_estimate);

typedef SLEQP_RETCODE (*SLEQP_FACTORIZATION_FREE)(void** star);

typedef struct
{
  SLEQP_FACTORIZATION_SET_MATRIX set_matrix;
  SLEQP_FACTORIZATION_SOLVE solve;
  SLEQP_FACTORIZATION_SOLUTION solution;
  SLEQP_FACTORIZATION_CONDITION_ESTIMATE condition_estimate;
  SLEQP_FACTORIZATION_FREE free;
} SleqpFactorizationCallbacks;

#endif /* SLEQP_FACTORIZATION
_TYPES_H */
