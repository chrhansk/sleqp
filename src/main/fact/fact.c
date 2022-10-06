#include "fact.h"

#include <string.h>

#include "fail.h"
#include "log.h"
#include "mem.h"

struct SleqpFact
{
  int refcount;

  char* name;
  char* version;

  SleqpFactCallbacks callbacks;
  SLEQP_FACT_FLAGS flags;
  void* fact_data;
};

SLEQP_RETCODE
sleqp_fact_create(SleqpFact** star,
                  const char* name,
                  const char* version,
                  SleqpParams* params,
                  SleqpFactCallbacks* callbacks,
                  SLEQP_FACT_FLAGS flags,
                  void* fact_data)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpFact* factorization = *star;

  *factorization = (SleqpFact){0};

  factorization->name    = strdup(name);
  factorization->version = strdup(version);

  factorization->refcount  = 1;
  factorization->callbacks = *callbacks;
  factorization->flags     = flags;
  factorization->fact_data = fact_data;

  return SLEQP_OKAY;
}

const char*
sleqp_fact_name(SleqpFact* factorization)
{
  return factorization->name;
}

const char*
sleqp_factorization_version(SleqpFact* factorization)
{
  return factorization->version;
}

SLEQP_RETCODE
sleqp_fact_set_matrix(SleqpFact* factorization, SleqpSparseMatrix* matrix)
{
#if SLEQP_DEBUG
  {
    if (factorization->flags & SLEQP_FACT_FLAGS_LOWER)
    {
      assert(sleqp_sparse_matrix_is_lower(matrix));
    }
  }
#endif

  SLEQP_CALL(
    factorization->callbacks.set_matrix(factorization->fact_data, matrix));

  return SLEQP_OKAY;
}

SLEQP_FACT_FLAGS
sleqp_fact_flags(SleqpFact* factorization)
{
  return factorization->flags;
}

SLEQP_RETCODE
sleqp_fact_solve(SleqpFact* factorization, const SleqpVec* rhs)
{
  SLEQP_CALL(factorization->callbacks.solve(factorization->fact_data, rhs));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_fact_solution(SleqpFact* factorization,
                    SleqpVec* sol,
                    int begin,
                    int end,
                    double zero_eps)
{
  SLEQP_CALL(factorization->callbacks
               .solution(factorization->fact_data, sol, begin, end, zero_eps));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_fact_cond(SleqpFact* factorization, double* condition)
{
  if (factorization->callbacks.condition)
  {
    SLEQP_CALL(
      factorization->callbacks.condition(factorization->fact_data, condition));
  }
  else
  {
    *condition = SLEQP_NONE;
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_fact_capture(SleqpFact* factorization)
{
  ++(factorization->refcount);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
factorization_free(SleqpFact** star)
{
  SleqpFact* factorization = *star;

  SLEQP_CALL(factorization->callbacks.free(&(factorization->fact_data)));

  sleqp_free(&factorization->version);
  sleqp_free(&factorization->name);

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_fact_release(SleqpFact** star)
{
  SleqpFact* factorization = *star;

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
