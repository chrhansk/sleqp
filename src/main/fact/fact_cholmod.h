#ifndef SLEQP_FACT_CHOLMOD_H
#define SLEQP_FACT_CHOLMOD_H

/**
 * @file fact_umfpack.h
 * @brief Defintion of CHOLMOD sparse factorization method.
 **/

#include "fact.h"
#include "types.h"

#include "sparse/sparse_matrix.h"

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_fact_cholmod_create(SleqpFact** star, SleqpParams* params);

#endif /* SLEQP_FACT_CHOLMOD_H */
