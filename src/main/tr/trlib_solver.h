#ifndef SLEQP_TRLIB_SOLVER_H
#define SLEQP_TRLIB_SOLVER_H

#include "tr_solver.h"

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_trlib_solver_create(SleqpTRSolver** star,
                          SleqpProblem* problem,
                          SleqpSettings* settings);

#endif /* SLEQP_TRLIB_SOLVER_H */
