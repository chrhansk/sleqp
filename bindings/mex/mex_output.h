#ifndef SLEQP_MEX_OUTPUT_H
#define SLEQP_MEX_OUTPUT_H

#include <mex.h>

#include "sleqp.h"

SLEQP_RETCODE
create_mex_output(SleqpProblem* problem,
                  SleqpSolver* solver,
                  mxArray** sol_star,
                  mxArray** info_star);

#endif /* SLEQP_MEX_OUTPUT_H */
