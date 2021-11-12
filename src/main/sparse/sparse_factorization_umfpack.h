#ifndef SLEQP_SPARSE_FACTORIZATION_UMFPACK_H
#define SLEQP_SPARSE_FACTORIZATION_UMFPACK_H

/**
 * @file sparse_factorization_umfpack.h
 * @brief Defintion of UMFPACK sparse factorization method.
 **/

#include "sparse_factorization.h"
#include "sparse_matrix.h"
#include "types.h"

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_sparse_factorization_umfpack_create(SleqpSparseFactorization** star,
                                          SleqpParams* params);

#endif /* SLEQP_SPARSE_FACTORIZATION_UMFPACK_H */
