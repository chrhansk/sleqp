#ifndef SLEQP_MEX_LSQ_H
#define SLEQP_MEX_LSQ_H

#include <mex.h>

#include "sleqp.h"

SLEQP_RETCODE
mex_lsq_func_create(SleqpFunc** star,
                    const mxArray* mex_x0,
                    const mxArray* mex_callbacks,
                    SleqpParams* params,
                    int num_variables,
                    int num_constraints);

#endif /* SLEQP_MEX_LSQ_H */
