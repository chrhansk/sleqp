#ifndef SLEQP_LSQR_TYPES_H
#define SLEQP_LSQR_TYPES_H

#include "sparse/vec.h"

typedef SLEQP_RETCODE (*SLEQP_LSQR_PROD_FORWARD)(const SleqpVec* direction,
                                                 SleqpVec* product,
                                                 void* data);

typedef SLEQP_RETCODE (*SLEQP_LSQR_PROD_ADJOINT)(const SleqpVec* direction,
                                                 SleqpVec* product,
                                                 void* data);

typedef struct
{
  SLEQP_LSQR_PROD_FORWARD prod_forward;
  SLEQP_LSQR_PROD_ADJOINT prod_adjoint;
} SleqpLSQRCallbacks;

#endif /* SLEQP_LSQR_TYPES_H */
