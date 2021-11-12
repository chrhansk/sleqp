#include "sparse_factorization.h"

#include <string.h>

#include "log.h"
#include "mem.h"

struct SleqpSparseFactorization
{
  int refcount;

  char* name;
  char* version;

  SleqpSparseFactorizationCallbacks callbacks;
  void* factorization_data;
};

SLEQP_RETCODE
sleqp_sparse_factorization_create(SleqpSparseFactorization** star,
                                  const char* name,
                                  const char* version,
                                  SleqpParams* params,
                                  SleqpSparseFactorizationCallbacks* callbacks,
                                  void* factorization_data)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpSparseFactorization* sparse_factorization = *star;

  *sparse_factorization = (SleqpSparseFactorization){0};

  sparse_factorization->name    = strdup(name);
  sparse_factorization->version = strdup(version);

  sparse_factorization->refcount           = 1;
  sparse_factorization->callbacks          = *callbacks;
  sparse_factorization->factorization_data = factorization_data;

  return SLEQP_OKAY;
}

const char*
sleqp_sparse_factorization_get_name(
  SleqpSparseFactorization* sparse_factorization)
{
  return sparse_factorization->name;
}

const char*
sleqp_sparse_factorization_get_version(
  SleqpSparseFactorization* sparse_factorization)
{
  return sparse_factorization->version;
}

SLEQP_RETCODE
sleqp_sparse_factorization_set_matrix(
  SleqpSparseFactorization* sparse_factorization,
  SleqpSparseMatrix* matrix)
{
  SLEQP_CALL(sparse_factorization->callbacks.set_matrix(
    sparse_factorization->factorization_data,
    matrix));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_sparse_factorization_solve(SleqpSparseFactorization* sparse_factorization,
                                 SleqpSparseVec* rhs)
{
  SLEQP_CALL(sparse_factorization->callbacks.solve(
    sparse_factorization->factorization_data,
    rhs));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_sparse_factorization_get_sol(
  SleqpSparseFactorization* sparse_factorization,
  SleqpSparseVec* sol,
  int begin,
  int end,
  double zero_eps)
{
  SLEQP_CALL(sparse_factorization->callbacks.get_sol(
    sparse_factorization->factorization_data,
    sol,
    begin,
    end,
    zero_eps));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_sparse_factorization_get_condition_estimate(
  SleqpSparseFactorization* sparse_factorization,
  double* condition_estimate)
{
  SLEQP_CALL(sparse_factorization->callbacks.get_condition_estimate(
    sparse_factorization->factorization_data,
    condition_estimate));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_sparse_factorization_capture(
  SleqpSparseFactorization* sparse_factorization)
{
  ++(sparse_factorization->refcount);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
sparse_factorization_free(SleqpSparseFactorization** star)
{
  SleqpSparseFactorization* sparse_factorization = *star;

  SLEQP_CALL(sparse_factorization->callbacks.free(
    &(sparse_factorization->factorization_data)));

  sleqp_free(&sparse_factorization->version);
  sleqp_free(&sparse_factorization->name);

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_sparse_factorization_release(SleqpSparseFactorization** star)
{
  SleqpSparseFactorization* sparse_factorization = *star;

  if (!sparse_factorization)
  {
    return SLEQP_OKAY;
  }

  if (--(sparse_factorization->refcount) == 0)
  {
    SLEQP_CALL(sparse_factorization_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
