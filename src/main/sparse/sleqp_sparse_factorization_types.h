#ifndef SLEQP_SPARSE_FACTORIZATION_TYPES_H
#define SLEQP_SPARSE_FACTORIZATION_TYPES_H

#include "sparse/sleqp_sparse_vec.h"
#include "sparse/sleqp_sparse_matrix.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpSparseFactorization SleqpSparseFactorization;

  typedef SLEQP_RETCODE
  (*SLEQP_SPARSE_FACTORIZATION_SET_MATRIX)(void* factorization_data,
                                           SleqpSparseMatrix* matrix);

  typedef SLEQP_RETCODE
  (*SLEQP_SPARSE_FACTORIZATION_SOLVE)(void* factorization_data,
                                      SleqpSparseVec* rhs);

  typedef SLEQP_RETCODE
  (*SLEQP_SPARSE_FACTORIZATION_GET_SOL)(void* factorization_data,
                                        SleqpSparseVec* sol,
                                        int begin,
                                        int end,
                                        double zero_eps);

  typedef SLEQP_RETCODE
  (*SLEQP_SPARSE_FACTORIZATION_GET_CONDITION_ESTIMATE)(void* factorization_data,
                                                       double* condition_estimate);

  typedef SLEQP_RETCODE
  (*SLEQP_SPARSE_FACTORIZATION_FREE)(void **star);

  typedef struct {
    SLEQP_SPARSE_FACTORIZATION_SET_MATRIX set_matrix;
    SLEQP_SPARSE_FACTORIZATION_SOLVE solve;
    SLEQP_SPARSE_FACTORIZATION_GET_SOL get_sol;
    SLEQP_SPARSE_FACTORIZATION_GET_CONDITION_ESTIMATE get_condition_estimate;
    SLEQP_SPARSE_FACTORIZATION_FREE free;
  } SleqpSparseFactorizationCallbacks;

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_SPARSE_FACTORIZATION_TYPES_H */
