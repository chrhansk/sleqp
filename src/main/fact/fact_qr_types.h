#ifndef SLEQP_FACT_QR_TYPES_H
#define SLEQP_FACT_QR_TYPES_H

#include "sparse/sparse_matrix.h"

typedef struct SleqpFactQR SleqpFactQR;

typedef SLEQP_RETCODE (*SLEQP_QR_SET_MATRIX)(void* fact_data,
                                             SleqpSparseMatrix* matrix);

typedef SLEQP_RETCODE (*SLEQP_QR_SOLVE_TRI)(void* fact_data,
                                            const SleqpVec* rhs,
                                            SleqpVec* sol);

typedef SLEQP_RETCODE (*SLEQP_QR_SOLVE_TRI_TRANS)(void* fact_data,
                                                  const SleqpVec* rhs,
                                                  SleqpVec* sol);

typedef SLEQP_RETCODE (*SLEQP_QR_MULT_ORTH)(void* fact_data,
                                            const SleqpVec* direction,
                                            SleqpVec* product);

typedef SLEQP_RETCODE (*SLEQP_QR_MULT_ORTH_TRANS)(void* fact_data,
                                                  const SleqpVec* direction,
                                                  SleqpVec* product);

typedef SLEQP_RETCODE (*SLEQP_QR_FREE)(void** star);

typedef struct
{
  SLEQP_QR_SET_MATRIX set_matrix;
  SLEQP_QR_SOLVE_TRI solve_tri;
  SLEQP_QR_SOLVE_TRI_TRANS solve_tri_trans;
  SLEQP_QR_MULT_ORTH mult_orth;
  SLEQP_QR_MULT_ORTH_TRANS mult_orth_trans;
  SLEQP_QR_FREE free;
} SleqpQRCallbacks;

#endif /* SLEQP_FACT_QR_TYPES_H */
