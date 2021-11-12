#ifndef SLEQP_LSQR_TYPES_H
#define SLEQP_LSQR_TYPES_H

#include "sparse/sparse_vec.h"

typedef SLEQP_RETCODE (*SLEQP_LSQR_PROD_FORWARD)(const SleqpSparseVec* direction,
                                                 SleqpSparseVec* product,
                                                 void* data);

typedef SLEQP_RETCODE (*SLEQP_LSQR_PROD_ADJOINT)(const SleqpSparseVec* direction,
                                                 SleqpSparseVec* product,
                                                 void* data);

typedef struct {
        SLEQP_LSQR_PROD_FORWARD prod_forward;
        SLEQP_LSQR_PROD_ADJOINT prod_adjoint;
} SleqpLSQRCallbacks;

#endif /* SLEQP_LSQR_TYPES_H */
