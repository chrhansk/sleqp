#ifndef SLEQP_FACT_MUMPS_H
#define SLEQP_FACT_MUMPS_H

#include "fact.h"
#include "sparse/mat.h"
#include "types.h"

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_fact_mumps_create(SleqpFact** star, SleqpSettings* settings);

#endif /* SLEQP_FACT_MUMPS_H */
