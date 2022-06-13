#ifndef SLEQP_FACT_MUMPS_H
#define SLEQP_FACT_MUMPS_H

#include "fact.h"
#include "sparse/sparse_matrix.h"
#include "types.h"

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_fact_mumps_create(SleqpFact** star, SleqpParams* params);

#endif /* SLEQP_FACT_MUMPS_H */
