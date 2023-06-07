#ifndef SLEQP_FIXED_VAR_FUNC_H
#define SLEQP_FIXED_VAR_FUNC_H

/**
 * @file fixed_var_func.h
 * @brief A function adapting to a set of fixed variables
 **/

#include "func.h"
#include "settings.h"

SLEQP_RETCODE
sleqp_fixed_var_func_create(SleqpFunc** star,
                            SleqpFunc* func,
                            int num_fixed,
                            const int* fixed_indices,
                            const double* fixed_values);

SLEQP_RETCODE
sleqp_fixed_var_lsq_func_create(SleqpFunc** star,
                                SleqpFunc* func,
                                SleqpSettings* settings,
                                int num_fixed,
                                const int* fixed_indices,
                                const double* fixed_values);

SLEQP_RETCODE
sleqp_fixed_var_dyn_func_create(SleqpFunc** star,
                                SleqpFunc* func,
                                int num_fixed,
                                const int* fixed_indices,
                                const double* fixed_values);

#endif /* SLEQP_FIXED_VAR_FUNC_H */
