#ifndef MEX_FUNC_H
#define MEX_FUNC_H

#include <mex.h>

#include "sleqp.h"

SLEQP_RETCODE
mex_func_create(SleqpFunc** star,
                const mxArray* ptr,
                SleqpParams* params,
                int num_variables,
                int num_constraints);

#endif /* MEX_FUNC_H */
