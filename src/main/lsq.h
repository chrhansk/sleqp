#ifndef SLEQP_LSQ_H
#define SLEQP_LSQ_H

#include "pub_lsq.h"

#ifdef __cplusplus
extern "C" {
#endif


  SLEQP_EXPORT
  double sleqp_lsq_func_get_levenberg_marquardt(SleqpFunc* func);

  SLEQP_EXPORT
  int sleqp_lsq_func_num_residuals(SleqpFunc* func);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_lsq_func_residuals(SleqpFunc* func,
                                         SleqpSparseVec* residuals);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_lsq_func_jac_forward(SleqpFunc* func,
                                           const SleqpSparseVec* forward_direction,
                                           SleqpSparseVec* product);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_lsq_func_jac_adjoint(SleqpFunc* func,
                                           const SleqpSparseVec* adjoint_direction,
                                           SleqpSparseVec* product);

#ifdef __cplusplus
}
#endif


#endif /* SLEQP_LSQ_H */
