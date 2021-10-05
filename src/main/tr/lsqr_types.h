#ifndef SLEQP_LSQR_TYPES_H
#define SLEQP_LSQR_TYPES_H

#include "sparse/sparse_vec.h"

#ifdef __cplusplus
extern "C" {
#endif

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

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_LSQR_TYPES_H */
