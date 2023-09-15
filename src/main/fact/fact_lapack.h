#ifndef SLEQP_FACT_LAPACK_H
#define SLEQP_FACT_LAPACK_H

#include "fact.h"
#include "sparse/mat.h"
#include "types.h"

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_fact_lapack_create(SleqpFact** star, SleqpSettings* settings);

#endif /* SLEQP_FACT_LAPACK_H */
