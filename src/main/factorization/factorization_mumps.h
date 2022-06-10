#ifndef SLEQP_FACTORIZATION_MUMPS_H
#define SLEQP_FACTORIZATION_MUMPS_H

#include "factorization.h"
#include "sparse/sparse_matrix.h"
#include "types.h"

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_fact_mumps_create(SleqpFact** star, SleqpParams* params);

#endif /* SLEQP_FACTORIZATION_MUMPS_H */
