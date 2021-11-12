#ifndef SLEQP_SPARSE_FACTORIZATION_MUMPS_H
#define SLEQP_SPARSE_FACTORIZATION_MUMPS_H

#include "sparse_factorization.h"
#include "sparse_matrix.h"
#include "types.h"

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_sparse_factorization_mumps_create(SleqpSparseFactorization** star,
                                        SleqpParams* params);

#endif /* SLEQP_SPARSE_FACTORIZATION_MUMPS_H */
