#ifndef SLEQP_SPARSE_FACTORIZATION_H
#define SLEQP_SPARSE_FACTORIZATION_H

#include "sleqp_types.h"

#include "sleqp_sparse_matrix.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpSparseFactorization SleqpSparseFactorization;

  SLEQP_RETCODE sleqp_sparse_factorization_create(SleqpSparseFactorization** star,
                                                  SleqpSparseMatrix* matrix);

  SLEQP_RETCODE sleqp_sparse_factorization_solve(SleqpSparseFactorization* factorization,
                                                 SleqpSparseVec* rhs);

  SLEQP_RETCODE sleqp_sparse_factorization_get_sol(SleqpSparseFactorization* factorization,
                                                   SleqpSparseVec* sol,
                                                   int begin,
                                                   int end);

  SLEQP_RETCODE sleqp_sparse_factorization_free(SleqpSparseFactorization** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_SPARSE_FACTORIZATION_H */
