#ifndef SLEQP_CUTEST_CONSTRAINED_H
#define SLEQP_CUTEST_CONSTRAINED_H

#include "sleqp.h"

#include "sleqp_cutest_data.h"
#include "sleqp_cutest_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_cutest_cons_func_create(SleqpFunc** star,
                                              int num_variables,
                                              int num_constraints,
                                              int num_linear,
                                              SleqpParams* params);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_cutest_cons_problem_create(SleqpProblem** star,
                                                 SleqpCutestData* data,
                                                 SleqpParams* params,
                                                 bool force_nonlinear);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_cutest_linear_offset(SleqpFunc* func,
                                           SleqpSparseVec* offset);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_cutest_eval_linear(SleqpFunc* func,
                                         SleqpSparseMatrix* coeffs);


#ifdef __cplusplus
}
#endif

#endif /* SLEQP_CUTEST_CONSTRAINED_H */
