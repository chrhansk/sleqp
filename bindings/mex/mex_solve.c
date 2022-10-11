#include "mex_solve.h"

#include <assert.h>

#include "mex_solve_common.h"

SLEQP_RETCODE
mex_command_solve(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  assert(nlhs == 2);
  assert(nrhs == 3);

  const mxArray* mex_x0      = prhs[0];
  const mxArray* mex_funcs   = prhs[1];
  const mxArray* mex_options = prhs[2];

  SLEQP_CALL(mex_solve(&plhs[0],
                       &plhs[1],
                       SLEQP_FUNC_TYPE_REGULAR,
                       mex_x0,
                       mex_funcs,
                       mex_options));

  return SLEQP_OKAY;
}
