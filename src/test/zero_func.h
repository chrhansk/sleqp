#ifndef ZERO_FUNC_H
#define ZERO_FUNC_H

#include "func.h"
#include "settings.h"

SLEQP_RETCODE
zero_func_create(SleqpFunc** star, int num_variables, int num_constraints);

SLEQP_RETCODE
zero_lsq_func_create(SleqpFunc** star,
                     SleqpSettings* settings,
                     int num_variables,
                     int num_constraints,
                     int num_residuals);

#endif /* ZERO_FUNC_H */
