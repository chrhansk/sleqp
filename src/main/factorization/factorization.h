#ifndef SLEQP_FACTORIZATION_H
#define SLEQP_FACTORIZATION_H

#include "factorization_types.h"
#include "params.h"
#include "timer.h"
#include "types.h"

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_factorization_create(SleqpFactorization** star,
                           const char* name,
                           const char* version,
                           SleqpParams* params,
                           SleqpFactorizationCallbacks* callbacks,
                           void* factorization_data);

const char*
sleqp_factorization_name(SleqpFactorization* factorization);

const char*
sleqp_factorization_version(SleqpFactorization* factorization);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_factorization_create_default(SleqpFactorization** star,
                                   SleqpParams* params);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_factorization_capture(SleqpFactorization* factorization);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_factorization_set_matrix(SleqpFactorization* factorization,
                               SleqpSparseMatrix* matrix);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_factorization_solve(SleqpFactorization* factorization,
                          SleqpSparseVec* rhs);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_factorization_solution(SleqpFactorization* factorization,
                             SleqpSparseVec* sol,
                             int begin,
                             int end,
                             double zero_eps);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_factorization_condition_estimate(SleqpFactorization* factorization,
                                       double* condition_estimate);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_factorization_capture(SleqpFactorization* factorization);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_factorization_release(SleqpFactorization** star);

#endif /* SLEQP_FACTORIZATION_H */
