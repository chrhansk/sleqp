#ifndef SLEQP_MEX_PROBLEM_H
#define SLEQP_MEX_PROBLEM_H

#include <mex.h>

#include "sleqp.h"

SLEQP_RETCODE
mex_problem_create(SleqpProblem** star,
                   SleqpSettings* settings,
                   SLEQP_FUNC_TYPE func_type,
                   const mxArray* mex_x0,
                   const mxArray* mex_funcs,
                   const mxArray* mex_options);

SLEQP_RETCODE
mex_create_vec_from_array(SleqpVec** star, const mxArray* array);

#endif /* SLEQP_MEX_PROBLEM_H */
