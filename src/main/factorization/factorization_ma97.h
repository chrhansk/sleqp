#ifndef SLEQP_FACTORIZATION_MA97_H
#define SLEQP_FACTORIZATION_MA97_H

#include "factorization.h"
#include "sparse/sparse_matrix.h"
#include "types.h"

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_fact_ma97_create(SleqpFact** star, SleqpParams* params);

#endif /* SLEQP_FACTORIZATION_MA97_H */
