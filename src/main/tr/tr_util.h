#ifndef SLEQP_TR_UTIL_H
#define SLEQP_TR_UTIL_H

#include "sparse/vec.h"
#include "settings.h"

SLEQP_RETCODE
sleqp_tr_compute_bdry_sol(const SleqpVec* previous,
                          const SleqpVec* direction,
                          SleqpSettings* settings,
                          double radius,
                          SleqpVec* result);

#endif /* SLEQP_TR_UTIL_H */
