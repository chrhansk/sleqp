#ifndef SLEQP_MEX_FUNC_H
#define SLEQP_MEX_FUNC_H

#include <mex.h>

#include "sleqp.h"

SLEQP_RETCODE
mex_func_create(SleqpFunc** star,
                const mxArray* mex_callbacks,
                SleqpSettings* settings,
                int num_vars,
                int num_cons,
                bool with_hess);

#endif /* SLEQP_MEX_FUNC_H */
