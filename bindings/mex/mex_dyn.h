#ifndef SLEQP_MEX_DYN_H
#define SLEQP_MEX_DYN_H

#include <mex.h>

#include "sleqp.h"

SLEQP_RETCODE
mex_dyn_func_create(SleqpFunc** star,
                    const mxArray* mex_x0,
                    const mxArray* mex_callbacks,
                    SleqpSettings* settings,
                    int num_vars,
                    int num_cons,
                    bool with_hess);

#endif /* SLEQP_MEX_DYN_H */
