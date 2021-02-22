#ifndef SLEQP_SPARSE_FACTORIZATION_UMFPACK_H
#define SLEQP_SPARSE_FACTORIZATION_UMFPACK_H

/**
 * @file sleqp_sparse_factorization_umfpack.h
 * @brief Defintion of UMFPACK sparse factorization method.
 **/

#include "sleqp_types.h"
#include "sleqp_sparse_matrix.h"
#include "sleqp_sparse_factorization.h"

#ifdef __cplusplus
extern "C" {
#endif

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_sparse_factorization_umfpack_create(SleqpSparseFactorization** star,
                                                          SleqpParams* params);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_SPARSE_FACTORIZATION_UMFPACK_H */
