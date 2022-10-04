#ifndef SLEQP_MEX_HESS_H
#define SLEQP_MEX_HESS_H

#include <mex.h>

#include "sleqp.h"

typedef struct
{
  struct
  {
    mxArray* hess;
  } callbacks;

  SleqpParams* params;

  mxArray* hess_dir;

  double* direction;
  double* product;

  mxArray* obj_dual;
  mxArray* cons_dual;

  mxArray** args;
  int nargs;

} MexHess;

SLEQP_RETCODE
mex_hess_init(MexHess* hess,
              SleqpParams* params,
              const mxArray* mex_callbacks,
              const int num_vars,
              const int num_cons);

SLEQP_RETCODE
mex_hess_prod(MexHess* hess,
              mxArray* primal,
              const double* obj_dual,
              const SleqpVec* direction,
              const SleqpVec* cons_duals,
              mxArray** rhs,
              int nrhs,
              SleqpVec* result);

SLEQP_RETCODE
mex_hess_free(MexHess* hess);

#endif /* SLEQP_MEX_HESS_H */
