#ifndef SLEQP_CUTEST_UNCONSTRAINED_H
#define SLEQP_CUTEST_UNCONSTRAINED_H

#include "sleqp.h"

#include "sleqp_cutest_data.h"

#ifdef __cplusplus
extern "C" {
#endif

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_cutest_uncons_func_create(SleqpFunc** star,
                                                int num_variables,
                                                SleqpParams* params);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_cutest_uncons_problem_create(SleqpProblem** star,
                                                   SleqpCutestData* data,
                                                   SleqpParams* params);


#ifdef __cplusplus
}
#endif

#endif /* SLEQP_CUTEST_UNCONSTRAINED_H */
