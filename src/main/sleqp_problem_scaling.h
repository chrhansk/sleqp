#ifndef SLEQP_PROBLEM_SCALING_H
#define SLEQP_PROBLEM_SCALING_H

#include "sleqp_scale.h"
#include "sleqp_options.h"
#include "sleqp_params.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpProblemScaling SleqpProblemScaling;

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_problem_scaling_create(SleqpProblemScaling** problem_scaling,
                                             SleqpScaling* scaling_data,
                                             SleqpProblem* problem,
                                             SleqpParams* params,
                                             SleqpOptions* options);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_problem_scaling_flush(SleqpProblemScaling* problem_scaling);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_problem_scaling_capture(SleqpProblemScaling* scaling);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_problem_scaling_release(SleqpProblemScaling** star);

  SleqpProblem* sleqp_problem_scaling_get_problem(SleqpProblemScaling* scaling);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_PROBLEM_SCALING_H */
