#ifndef SLEQP_STEIHAUG_SOLVER_H
#define SLEQP_STEIHAUG_SOLVER_H

#include "sleqp_tr_solver.h"

#ifdef __cplusplus
extern "C" {
#endif

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_steihaug_solver_create(SleqpTRSolver** star,
                                             SleqpProblem* problem,
                                             SleqpParams* params,
                                             SleqpOptions* options);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_STEIHAUG_SOLVER_H */
