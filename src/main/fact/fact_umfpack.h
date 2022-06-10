#ifndef SLEQP_FACT_UMFPACK_H
#define SLEQP_FACT_UMFPACK_H

/**
 * @file fact_umfpack.h
 * @brief Defintion of UMFPACK sparse factorization method.
 **/

#include "fact/fact.h"
#include "types.h"

#include "sparse/sparse_matrix.h"

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_fact_umfpack_create(SleqpFact** star, SleqpParams* params);

#endif /* SLEQP_FACT_UMFPACK_H */
