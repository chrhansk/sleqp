#ifndef SLEQP_FACTORIZATION_MA86_H
#define SLEQP_FACTORIZATION_MA86_H

#include "factorization.h"
#include "sparse/sparse_matrix.h"
#include "types.h"

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_factorization_ma86_create(SleqpFactorization** star, SleqpParams* params);

#endif /* SLEQP_FACTORIZATION_MA86_H */
