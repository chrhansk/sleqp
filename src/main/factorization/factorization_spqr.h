#ifndef SLEQP_FACTORIZATION_SPQR_H
#define SLEQP_FACTORIZATION_SPQR_H

/**
 * @file factorization_umfpack.h
 * @brief Defintion of SPQR sparse factorization method.
 **/

#include "factorization.h"
#include "types.h"

#include "sparse/sparse_matrix.h"

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_fact_spqr_create(SleqpFact** star, SleqpParams* params);

#endif /* SLEQP_FACTORIZATION_SPQR_H */
