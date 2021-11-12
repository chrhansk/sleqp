#ifndef SLEQP_TR_UTIL_H
#define SLEQP_TR_UTIL_H

#include "params.h"
#include "sparse/sparse_vec.h"

SLEQP_RETCODE sleqp_tr_compute_bdry_sol(const SleqpSparseVec* previous,
                                        const SleqpSparseVec* direction,
                                        SleqpParams* params,
                                        double radius,
                                        SleqpSparseVec* result);

#endif /* SLEQP_TR_UTIL_H */
