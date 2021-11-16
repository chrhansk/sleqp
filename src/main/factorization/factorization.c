#include "factorization.h"

#include <string.h>

#include "log.h"
#include "mem.h"

struct SleqpFactorization
{
  int refcount;

  char* name;
  char* version;

  SleqpFactorizationCallbacks callbacks;
  void* factorization_data;
};

SLEQP_RETCODE
sleqp_factorization_create(SleqpFactorization** star,
                           const char* name,
                           const char* version,
                           SleqpParams* params,
                           SleqpFactorizationCallbacks* callbacks,
                           void* factorization_data)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpFactorization* factorization = *star;

  *factorization = (SleqpFactorization){0};

  factorization->name    = strdup(name);
  factorization->version = strdup(version);

  factorization->refcount           = 1;
  factorization->callbacks          = *callbacks;
  factorization->factorization_data = factorization_data;

  return SLEQP_OKAY;
}

const char*
sleqp_factorization_name(SleqpFactorization* factorization)
{
  return factorization->name;
}

const char*
sleqp_factorization_version(SleqpFactorization* factorization)
{
  return factorization->version;
}

SLEQP_RETCODE
sleqp_factorization_set_matrix(SleqpFactorization* factorization,
                               SleqpSparseMatrix* matrix)
{
  SLEQP_CALL(
    factorization->callbacks.set_matrix(factorization->factorization_data,
                                        matrix));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_factorization_solve(SleqpFactorization* factorization,
                          SleqpSparseVec* rhs)
{
  SLEQP_CALL(
    factorization->callbacks.solve(factorization->factorization_data, rhs));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_factorization_solution(SleqpFactorization* factorization,
                             SleqpSparseVec* sol,
                             int begin,
                             int end,
                             double zero_eps)
{
  SLEQP_CALL(
    factorization->callbacks
      .solution(factorization->factorization_data, sol, begin, end, zero_eps));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_factorization_condition_estimate(SleqpFactorization* factorization,
                                       double* condition_estimate)
{
  SLEQP_CALL(factorization->callbacks.condition_estimate(
    factorization->factorization_data,
    condition_estimate));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_factorization_capture(SleqpFactorization* factorization)
{
  ++(factorization->refcount);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
factorization_free(SleqpFactorization** star)
{
  SleqpFactorization* factorization = *star;

  SLEQP_CALL(
    factorization->callbacks.free(&(factorization->factorization_data)));

  sleqp_free(&factorization->version);
  sleqp_free(&factorization->name);

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_factorization_release(SleqpFactorization** star)
{
  SleqpFactorization* factorization = *star;

  if (!factorization)
  {
    return SLEQP_OKAY;
  }

  if (--(factorization->refcount) == 0)
  {
    SLEQP_CALL(factorization_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
