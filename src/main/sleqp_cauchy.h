#ifndef SLEQP_CAUCHY_H
#define SLEQP_CAUCHY_H

#include "sleqp.h"
#include "lp/sleqp_lpi.h"

SLEQP_RETCODE sleqp_cauchy_direction(SleqpProblem* problem,
                                     SleqpSparseVec* x,
                                     SleqpLPi* lp_interface,
                                     SleqpSparseVec* func_grad,
                                     SleqpSparseMatrix* cons_jac);

#endif /* SLEQP_CAUCHY_H */
