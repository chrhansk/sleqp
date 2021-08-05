#ifndef SLEQP_FIXED_VAR_FUNC_H
#define SLEQP_FIXED_VAR_FUNC_H

/**
 * @file fixed_var_func.h
 * @brief A function adapting to a set of fixed variables
 **/

#include "func.h"
#include "params.h"

#ifdef __cplusplus
extern "C" {
#endif

  SLEQP_RETCODE sleqp_fixed_var_func_create(SleqpFunc** star,
                                            SleqpFunc* func,
                                            int num_fixed,
                                            const int* fixed_indices,
                                            const double* fixed_values);

  SLEQP_RETCODE sleqp_fixed_var_lsq_func_create(SleqpFunc** star,
                                                SleqpFunc* func,
                                                SleqpParams* params,
                                                int num_fixed,
                                                const int* fixed_indices,
                                                const double* fixed_values);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_FIXED_VAR_FUNC_H */
