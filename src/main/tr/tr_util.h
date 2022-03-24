#ifndef SLEQP_TR_UTIL_H
#define SLEQP_TR_UTIL_H

#include "params.h"
#include "sparse/vec.h"

SLEQP_RETCODE
sleqp_tr_compute_bdry_sol(const SleqpVec* previous,
                          const SleqpVec* direction,
                          SleqpParams* params,
                          double radius,
                          SleqpVec* result);

#endif /* SLEQP_TR_UTIL_H */
