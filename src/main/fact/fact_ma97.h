#ifndef SLEQP_FACT_MA97_H
#define SLEQP_FACT_MA97_H

#include "fact.h"
#include "sparse/mat.h"
#include "types.h"

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_fact_ma97_create(SleqpFact** star, SleqpParams* params);

#endif /* SLEQP_FACT_MA97_H */
