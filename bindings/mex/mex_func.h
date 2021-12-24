#ifndef SLEQP_MEX_FUNC_H
#define SLEQP_MEX_FUNC_H

#include <mex.h>

#include "sleqp.h"

SLEQP_RETCODE
mex_func_create(SleqpFunc** star,
                const mxArray* ptr,
                SleqpParams* params,
                int num_variables,
                int num_constraints);

#endif /* SLEQP_MEX_FUNC_H */
