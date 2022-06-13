#ifndef SLEQP_FACT_SPQR_H
#define SLEQP_FACT_SPQR_H

/**
 * @file fact_spqr.h
 * @brief Defintion of SPQR sparse factorization method.
 **/

#include "fact.h"
#include "types.h"

#include "sparse/sparse_matrix.h"

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_fact_spqr_create(SleqpFact** star, SleqpParams* params);

#endif /* SLEQP_FACT_SPQR_H */
