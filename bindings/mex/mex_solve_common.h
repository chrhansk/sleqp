#ifndef SLEQP_MEX_SOLVE_COMMON_H
#define SLEQP_MEX_SOLVE_COMMON_H

#include <mex.h>

#include "sleqp.h"

SLEQP_RETCODE
mex_solve(mxArray** sol_star,
          mxArray** info_star,
          bool lsq,
          const mxArray* mex_x0,
          const mxArray* mex_funcs,
          const mxArray* mex_options);

#endif /* SLEQP_MEX_SOLVE_COMMON_H */
