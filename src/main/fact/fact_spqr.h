#ifndef SLEQP_FACT_SPQR_H
#define SLEQP_FACT_SPQR_H

/**
 * @file fact_spqr.h
 * @brief Defintion of SPQR sparse factorization method.
 **/

#include "fact.h"
#include "fact_qr.h"

#include "sparse/mat.h"

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_fact_spqr_create(SleqpFactQR** star, SleqpSettings* settings);

#endif /* SLEQP_FACT_SPQR_H */
