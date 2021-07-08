#ifndef SLEQP_SPARSE_FACTORIZATION_MA97_H
#define SLEQP_SPARSE_FACTORIZATION_MA97_H

#include "types.h"
#include "sparse_matrix.h"
#include "sparse_factorization.h"

#ifdef __cplusplus
extern "C" {
#endif

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_sparse_factorization_ma97_create(SleqpSparseFactorization** star,
                                                       SleqpParams* params);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_SPARSE_FACTORIZATION_MA97_H */
