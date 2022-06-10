#ifndef SLEQP_FACT_H
#define SLEQP_FACT_H

#include "fact_types.h"
#include "params.h"
#include "timer.h"
#include "types.h"

typedef enum
{
  SLEQP_FACT_FLAGS_NONE  = 0,        /** Nothing **/
  SLEQP_FACT_FLAGS_PSD   = (1 << 0), /** Requires positive definiteness **/
  SLEQP_FACT_FLAGS_LOWER = (1 << 1)  /** Pass only lower triangular part **/
} SLEQP_FACT_FLAGS;

/**
 * Creates a new factorization. Factorizations are used to solve the
 * *symmetric* (but possible indefinite) systems occuring during the
 * optimization
 **/
SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_fact_create(SleqpFact** star,
                  const char* name,
                  const char* version,
                  SleqpParams* params,
                  SleqpFactorizationCallbacks* callbacks,
                  SLEQP_FACT_FLAGS flags,
                  void* fact_data);

const char*
sleqp_fact_name(SleqpFact* factorization);

const char*
sleqp_fact_version(SleqpFact* factorization);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_fact_set_matrix(SleqpFact* factorization, SleqpSparseMatrix* matrix);

SLEQP_FACT_FLAGS
sleqp_fact_flags(SleqpFact* factorization);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_fact_create_default(SleqpFact** star, SleqpParams* params);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_fact_capture(SleqpFact* factorization);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_fact_solve(SleqpFact* factorization, const SleqpVec* rhs);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_fact_solution(SleqpFact* factorization,
                    SleqpVec* sol,
                    int begin,
                    int end,
                    double zero_eps);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_fact_condition(SleqpFact* factorization, double* condition_estimate);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_fact_release(SleqpFact** star);

#endif /* SLEQP_FACT_H */
