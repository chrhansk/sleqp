#ifndef SLEQP_FACT_MA57_H
#define SLEQP_FACT_MA57_H

#include "fact.h"
#include "sparse/sparse_matrix.h"
#include "types.h"

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_fact_ma57_create(SleqpFact** star, SleqpParams* params);

#endif /* SLEQP_FACT_MA57_H */
