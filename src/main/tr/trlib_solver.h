#ifndef SLEQP_TRLIB_SOLVER_H
#define SLEQP_TRLIB_SOLVER_H

#include "tr_solver.h"

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_trlib_solver_create(SleqpTRSolver** star,
                          SleqpProblem* problem,
                          SleqpParams* params,
                          SleqpOptions* options);

#endif /* SLEQP_TRLIB_SOLVER_H */
