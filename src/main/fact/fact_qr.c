#include "fact_qr.h"

#include <assert.h>
#include <string.h>

#include "mem.h"

struct SleqpFactQR
{
  int refcount;

  char* name;
  char* version;

  int num_rows;
  int num_cols;

  SleqpQRCallbacks callbacks;

  void* fact_data;
};

SLEQP_RETCODE
sleqp_qr_create(SleqpFactQR** star,
                const char* name,
                const char* version,
                SleqpParams* params,
                SleqpQRCallbacks* callbacks,
                void* fact_data)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpFactQR* qr = *star;

  *qr = (SleqpFactQR){0};

  qr->name    = strdup(name);
  qr->version = strdup(version);

  qr->refcount  = 1;
  qr->callbacks = *callbacks;
  qr->fact_data = fact_data;

  return SLEQP_OKAY;
}

const char*
sleqp_qr_name(SleqpFactQR* fact)
{
  return fact->name;
}

const char*
sleqp_qr_version(SleqpFactQR* fact)
{
  return fact->version;
}

SLEQP_RETCODE
sleqp_qr_set_matrix(SleqpFactQR* fact, SleqpMat* matrix)
{
  fact->num_rows = sleqp_mat_num_rows(matrix);
  fact->num_cols = sleqp_mat_num_cols(matrix);

  assert(fact->num_rows >= fact->num_cols);

  SLEQP_CALL(fact->callbacks.set_matrix(fact->fact_data, matrix));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_qr_solve_tri(SleqpFactQR* fact, const SleqpVec* rhs, SleqpVec* sol)
{
  assert(rhs->dim == fact->num_cols);
  assert(sol->dim == fact->num_cols);

  SLEQP_CALL(fact->callbacks.solve_tri(fact->fact_data, rhs, sol));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_qr_solve_tri_trans(SleqpFactQR* fact, const SleqpVec* rhs, SleqpVec* sol)
{
  assert(rhs->dim == fact->num_cols);
  assert(sol->dim == fact->num_cols);

  SLEQP_CALL(fact->callbacks.solve_tri_trans(fact->fact_data, rhs, sol));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_qr_mult_orth(SleqpFactQR* fact,
                   const SleqpVec* direction,
                   SleqpVec* product)
{
  assert(direction->dim == fact->num_rows);
  assert(product->dim == fact->num_rows);

  SLEQP_CALL(fact->callbacks.mult_orth(fact->fact_data, direction, product));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_qr_mult_orth_trans(SleqpFactQR* fact,
                         const SleqpVec* direction,
                         SleqpVec* product)
{
  assert(product->dim == fact->num_rows);
  assert(direction->dim == fact->num_rows);

  SLEQP_CALL(
    fact->callbacks.mult_orth_trans(fact->fact_data, direction, product));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_qr_capture(SleqpFactQR* fact)
{
  ++(fact->refcount);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
qr_free(SleqpFactQR** star)
{
  SleqpFactQR* qr = *star;

  SLEQP_CALL(qr->callbacks.free(&qr->fact_data));

  sleqp_free(&(qr->version));
  sleqp_free(&(qr->name));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_qr_release(SleqpFactQR** star)
{
  SleqpFactQR* qr = *star;

  if (!qr)
  {
    return SLEQP_OKAY;
  }

  if (--(qr->refcount) == 0)
  {
    SLEQP_CALL(qr_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
