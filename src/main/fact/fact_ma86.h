#ifndef SLEQP_FACT_MA86_H
#define SLEQP_FACT_MA86_H

#include "fact.h"
#include "sparse/mat.h"
#include "types.h"

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_fact_ma86_create(SleqpFact** star, SleqpParams* params);

#endif /* SLEQP_FACT_MA86_H */
