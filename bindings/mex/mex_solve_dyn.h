#ifndef SLEQP_MEX_SOLVE_DYN_H
#define SLEQP_MEX_SOLVE_DYN_H

#include "sleqp.h"
#include <mex.h>

SLEQP_RETCODE
mex_command_solve_dyn(int nlhs,
                      mxArray* plhs[],
                      int nrhs,
                      const mxArray* prhs[]);

#endif /* SLEQP_MEX_SOLVE_DYN_H */
