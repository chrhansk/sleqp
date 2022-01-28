#include "lsq.h"

#include <assert.h>
#include <math.h>

#include "cmp.h"
#include "func.h"
#include "log.h"
#include "mem.h"
#include "sparse/sparse_matrix.h"

typedef struct SleqpLSQData
{
  int num_variables;
  int num_residuals;

  SleqpLSQCallbacks callbacks;

  double lm_factor;

  // num_residuals
  SleqpSparseVec* lsq_forward;
  SleqpSparseVec* lsq_residual;
  bool has_lsq_residual;

  // num_variables
  SleqpSparseVec* lsq_grad;
  SleqpSparseVec* lsq_hess_prod;

  double zero_eps;

  void* func_data;

} SleqpLSQData;

static SLEQP_RETCODE
lsq_func_set_value(SleqpFunc* func,
                   SleqpSparseVec* x,
                   SLEQP_VALUE_REASON reason,
                   bool* reject,
                   int* obj_grad_nnz,
                   int* cons_val_nnz,
                   int* cons_jac_nnz,
                   void* func_data)
{
  SleqpLSQData* lsq_data = (SleqpLSQData*)func_data;

  SLEQP_FUNC_CALL(lsq_data->callbacks.set_value(func,
                                                x,
                                                reason,
                                                reject,
                                                obj_grad_nnz,
                                                cons_val_nnz,
                                                cons_jac_nnz,
                                                lsq_data->func_data),
                  sleqp_func_has_flags(func, SLEQP_FUNC_INTERNAL),
                  "Error setting function value");

  lsq_data->has_lsq_residual = false;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_lsq_residual(SleqpFunc* func, SleqpLSQData* lsq_data)
{
  if (!(lsq_data->has_lsq_residual))
  {
    SLEQP_CALL(sleqp_sparse_vector_clear(lsq_data->lsq_residual));

    SLEQP_FUNC_CALL(lsq_data->callbacks.lsq_residuals(func,
                                                      lsq_data->lsq_residual,
                                                      lsq_data->func_data),
                    sleqp_func_has_flags(func, SLEQP_FUNC_INTERNAL),
                    "Error evaluating least squares residuals");

    lsq_data->has_lsq_residual = true;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
lsq_func_obj_val(SleqpFunc* func, double* obj_val, void* func_data)
{
  SleqpLSQData* lsq_data = (SleqpLSQData*)func_data;

  SLEQP_CALL(compute_lsq_residual(func, lsq_data));

  *obj_val = .5 * sleqp_sparse_vector_norm_sq(lsq_data->lsq_residual);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
lsq_func_obj_grad(SleqpFunc* func, SleqpSparseVec* obj_grad, void* func_data)
{
  SleqpLSQData* lsq_data = (SleqpLSQData*)func_data;

  SLEQP_CALL(sleqp_sparse_vector_clear(obj_grad));

  SLEQP_CALL(compute_lsq_residual(func, lsq_data));

  SLEQP_CALL(
    sleqp_lsq_func_jac_adjoint(func, lsq_data->lsq_residual, obj_grad));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
lsq_func_cons_val(SleqpFunc* func, SleqpSparseVec* cons_val, void* func_data)
{
  SleqpLSQData* lsq_data = (SleqpLSQData*)func_data;

  const int num_constraints = sleqp_func_num_cons(func);

  SLEQP_CALL(sleqp_sparse_vector_clear(cons_val));

  if (num_constraints != 0)
  {
    SLEQP_FUNC_CALL(
      lsq_data->callbacks.cons_val(func, cons_val, lsq_data->func_data),
      sleqp_func_has_flags(func, SLEQP_FUNC_INTERNAL),
      SLEQP_FUNC_ERROR_CONS_VAL);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
lsq_func_cons_jac(SleqpFunc* func, SleqpSparseMatrix* cons_jac, void* func_data)
{
  SleqpLSQData* lsq_data = (SleqpLSQData*)func_data;

  const int num_constraints = sleqp_func_num_cons(func);

  SLEQP_CALL(sleqp_sparse_matrix_clear(cons_jac));

  if (num_constraints != 0)
  {
    SLEQP_FUNC_CALL(
      lsq_data->callbacks.cons_jac(func, cons_jac, lsq_data->func_data),
      sleqp_func_has_flags(func, SLEQP_FUNC_INTERNAL),
      SLEQP_FUNC_ERROR_CONS_JAC);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_lsq_hess_prod(SleqpFunc* func,
                      SleqpLSQData* lsq_data,
                      const double* obj_dual,
                      const SleqpSparseVec* direction,
                      SleqpSparseVec* destination)
{
  if (obj_dual)
  {
    SLEQP_CALL(
      sleqp_lsq_func_jac_forward(func, direction, lsq_data->lsq_forward));

    SLEQP_CALL(
      sleqp_lsq_func_jac_adjoint(func, lsq_data->lsq_forward, destination));

    SLEQP_CALL(sleqp_sparse_vector_scale(destination, *obj_dual));
  }
  else
  {
    SLEQP_CALL(sleqp_sparse_vector_clear(destination));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
lsq_func_hess_product(SleqpFunc* func,
                      const double* obj_dual,
                      const SleqpSparseVec* direction,
                      const SleqpSparseVec* cons_duals,
                      SleqpSparseVec* product,
                      void* func_data)
{
  SleqpLSQData* lsq_data = (SleqpLSQData*)func_data;

  SLEQP_CALL(sleqp_sparse_vector_clear(lsq_data->lsq_forward));
  SLEQP_CALL(sleqp_sparse_vector_clear(lsq_data->lsq_hess_prod));

  const bool additional_term = (lsq_data->lm_factor != 0.);

  if (additional_term)
  {
    SleqpSparseVec* initial_product = lsq_data->lsq_hess_prod;

    SLEQP_CALL(compute_lsq_hess_prod(func,
                                     lsq_data,
                                     obj_dual,
                                     direction,
                                     initial_product));

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(initial_product,
                                              direction,
                                              1.,
                                              lsq_data->lm_factor,
                                              lsq_data->zero_eps,
                                              product));
  }
  else
  {
    SLEQP_CALL(
      compute_lsq_hess_prod(func, lsq_data, obj_dual, direction, product));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
lsq_func_free(void* func_data)
{
  if (!func_data)
  {
    return SLEQP_OKAY;
  }

  SleqpLSQData* lsq_data = (SleqpLSQData*)func_data;

  if (lsq_data->callbacks.func_free)
  {
    SLEQP_CALL(lsq_data->callbacks.func_free(lsq_data->func_data));
  }

  SLEQP_CALL(sleqp_sparse_vector_free(&lsq_data->lsq_hess_prod));

  SLEQP_CALL(sleqp_sparse_vector_free(&lsq_data->lsq_grad));

  SLEQP_CALL(sleqp_sparse_vector_free(&lsq_data->lsq_residual));

  SLEQP_CALL(sleqp_sparse_vector_free(&lsq_data->lsq_forward));

  sleqp_free(&lsq_data);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_lsq_func_create(SleqpFunc** fstar,
                      SleqpLSQCallbacks* callbacks,
                      int num_variables,
                      int num_constraints,
                      int num_residuals,
                      double lm_factor,
                      SleqpParams* params,
                      void* func_data)
{
  SleqpLSQData* data = NULL;

  SLEQP_CALL(sleqp_malloc(&data));

  *data = (SleqpLSQData){0};

  data->num_variables = num_variables;
  data->num_residuals = num_residuals;

  data->has_lsq_residual = false;

  data->callbacks = *callbacks;

  data->lm_factor = lm_factor;

  SLEQP_CALL(
    sleqp_sparse_vector_create_empty(&data->lsq_forward, num_residuals));

  SLEQP_CALL(
    sleqp_sparse_vector_create_empty(&data->lsq_residual, num_residuals));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->lsq_grad, num_variables));

  SLEQP_CALL(
    sleqp_sparse_vector_create_empty(&data->lsq_hess_prod, num_variables));

  data->zero_eps = sleqp_params_value(params, SLEQP_PARAM_ZERO_EPS);

  data->func_data = func_data;

  SleqpFuncCallbacks func_callbacks = {.set_value = lsq_func_set_value,
                                       .obj_val   = lsq_func_obj_val,
                                       .obj_grad  = lsq_func_obj_grad,
                                       .cons_val  = lsq_func_cons_val,
                                       .cons_jac  = lsq_func_cons_jac,
                                       .hess_prod = lsq_func_hess_product,
                                       .func_free = lsq_func_free};

  SLEQP_CALL(sleqp_func_create(fstar,
                               &func_callbacks,
                               num_variables,
                               num_constraints,
                               data));

  SleqpFunc* func = *fstar;

  SLEQP_FUNC_FLAGS flags = (SLEQP_FUNC_HESS_PSD | SLEQP_FUNC_HESS_INEXACT
                            | SLEQP_FUNC_HESS_INTERNAL);

  SLEQP_CALL(sleqp_func_flags_add(func, flags));

  SLEQP_CALL(sleqp_func_set_type(func, SLEQP_FUNC_TYPE_LSQ));

  return SLEQP_OKAY;
}

double
sleqp_lsq_func_get_levenberg_marquardt(SleqpFunc* func)
{
  assert(sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_LSQ);
  void* func_data = sleqp_func_get_data(func);
  assert(func_data);

  SleqpLSQData* lsq_data = (SleqpLSQData*)func_data;

  return lsq_data->lm_factor;
}

int
sleqp_lsq_func_num_residuals(SleqpFunc* func)
{
  assert(sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_LSQ);
  void* func_data = sleqp_func_get_data(func);
  assert(func_data);

  SleqpLSQData* lsq_data = (SleqpLSQData*)func_data;

  return lsq_data->num_residuals;
}

SLEQP_RETCODE
sleqp_lsq_func_residuals(SleqpFunc* func, SleqpSparseVec* residuals)
{
  assert(sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_LSQ);
  void* func_data = sleqp_func_get_data(func);
  assert(func_data);

  SleqpLSQData* lsq_data = (SleqpLSQData*)func_data;

  assert(residuals->dim == sleqp_lsq_func_num_residuals(func));

  SLEQP_CALL(compute_lsq_residual(func, lsq_data));

  SLEQP_CALL(sleqp_sparse_vector_copy(lsq_data->lsq_residual, residuals));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_lsq_func_jac_forward(SleqpFunc* func,
                           const SleqpSparseVec* forward_direction,
                           SleqpSparseVec* product)
{
  assert(sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_LSQ);
  void* func_data = sleqp_func_get_data(func);
  assert(func_data);

  SleqpLSQData* lsq_data = (SleqpLSQData*)func_data;

  assert(product != forward_direction);
  assert(forward_direction->dim == sleqp_func_num_vars(func));
  assert(product->dim == sleqp_lsq_func_num_residuals(func));

  SLEQP_CALL(sleqp_sparse_vector_clear(product));

  SLEQP_FUNC_CALL(lsq_data->callbacks.lsq_jac_forward(func,
                                                      forward_direction,
                                                      product,
                                                      lsq_data->func_data),
                  sleqp_func_has_flags(func, SLEQP_FUNC_INTERNAL),
                  "Error evaluating forward Jacobian product");

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_lsq_func_jac_adjoint(SleqpFunc* func,
                           const SleqpSparseVec* adjoint_direction,
                           SleqpSparseVec* product)
{
  assert(sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_LSQ);
  void* func_data = sleqp_func_get_data(func);
  assert(func_data);

  SleqpLSQData* lsq_data = (SleqpLSQData*)func_data;

  assert(product != adjoint_direction);
  assert(adjoint_direction->dim == sleqp_lsq_func_num_residuals(func));
  assert(product->dim == sleqp_func_num_vars(func));

  SLEQP_CALL(sleqp_sparse_vector_clear(product));

  SLEQP_FUNC_CALL(lsq_data->callbacks.lsq_jac_adjoint(func,
                                                      adjoint_direction,
                                                      product,
                                                      lsq_data->func_data),
                  sleqp_func_has_flags(func, SLEQP_FUNC_INTERNAL),
                  "Error evaluating adjoint Jacobian product");

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_lsq_func_set_callbacks(SleqpFunc* func, SleqpLSQCallbacks* callbacks)
{
  assert(sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_LSQ);
  void* func_data = sleqp_func_get_data(func);
  assert(func_data);

  SleqpLSQData* lsq_data = (SleqpLSQData*)func_data;

  lsq_data->callbacks = *callbacks;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_lsq_func_set_lm_factor(SleqpFunc* func, double lm_factor)
{
  assert(sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_LSQ);
  void* func_data = sleqp_func_get_data(func);
  assert(func_data);

  SleqpLSQData* lsq_data = (SleqpLSQData*)func_data;

  lsq_data->lm_factor = lm_factor;

  return SLEQP_OKAY;
}

void*
sleqp_lsq_func_get_data(SleqpFunc* func)
{
  assert(sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_LSQ);
  void* func_data = sleqp_func_get_data(func);
  assert(func_data);

  SleqpLSQData* lsq_data = (SleqpLSQData*)func_data;

  return lsq_data->func_data;
}
