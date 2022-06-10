#ifndef SLEQP_FACTORIZATION_H
#define SLEQP_FACTORIZATION_H

#include "factorization_types.h"
#include "params.h"
#include "timer.h"
#include "types.h"

typedef enum
{
  SLEQP_FACTORIZATION_NONE  = 0,        /** Nothing **/
  SLEQP_FACTORIZATION_PSD   = (1 << 0), /** Requires positive definiteness **/
  SLEQP_FACTORIZATION_LOWER = (1 << 1)  /** Pass only lower triangular part **/
} SLEQP_FACTORIZATION_FLAGS;

/**
 * Creates a new factorization. Factorizations are used to solve the
 * *symmetric* (but possible indefinite) systems occuring during the
 * optimization
 **/
SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_factorization_create(SleqpFactorization** star,
                           const char* name,
                           const char* version,
                           SleqpParams* params,
                           SleqpFactorizationCallbacks* callbacks,
                           SLEQP_FACTORIZATION_FLAGS flags,
                           void* factorization_data);

const char*
sleqp_factorization_name(SleqpFactorization* factorization);

const char*
sleqp_factorization_version(SleqpFactorization* factorization);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_factorization_set_matrix(SleqpFactorization* factorization,
                               SleqpSparseMatrix* matrix);

SLEQP_FACTORIZATION_FLAGS
sleqp_factorization_flags(SleqpFactorization* factorization);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_factorization_create_default(SleqpFactorization** star,
                                   SleqpParams* params);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_factorization_capture(SleqpFactorization* factorization);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_factorization_solve(SleqpFactorization* factorization,
                          const SleqpVec* rhs);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_factorization_solution(SleqpFactorization* factorization,
                             SleqpVec* sol,
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
