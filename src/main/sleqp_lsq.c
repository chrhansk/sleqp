#include "sleqp_lsq.h"

#include <assert.h>
#include <math.h>

#include "sleqp_cmp.h"
#include "sleqp_func.h"
#include "sleqp_mem.h"

typedef struct SleqpLSQData
{
  int num_variables;
  int num_residuals;

  SleqpLSQCallbacks callbacks;

  double levenberg_marquardt;

  // num_residuals
  SleqpSparseVec* lsq_forward;
  SleqpSparseVec* lsq_residual;
  bool has_lsq_residual;

  // num_variables
  SleqpSparseVec* lsq_hess_prod;
  SleqpSparseVec* lsq_grad;
  SleqpSparseVec* grad_cache;
  SleqpSparseVec* hess_prod;

  SleqpSparseVec* combined_hess_prod;

  double eps;

  void* func_data;

} SleqpLSQData;


static SLEQP_RETCODE lsq_func_set_value(SleqpFunc* func,
                                        SleqpSparseVec* x,
                                        SLEQP_VALUE_REASON reason,
                                        int* func_grad_nnz,
                                        int* cons_val_nnz,
                                        int* cons_jac_nnz,
                                        void* func_data)
{
  SleqpLSQData* lsq_data = (SleqpLSQData*) func_data;

  SLEQP_CALL(lsq_data->callbacks.set_value(func,
                                           x,
                                           reason,
                                           func_grad_nnz,
                                           cons_val_nnz,
                                           cons_jac_nnz,
                                           lsq_data->func_data));

  lsq_data->has_lsq_residual = false;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE compute_lsq_residual(SleqpFunc* func,
                                          SleqpLSQData* lsq_data)
{
  if(!(lsq_data->has_lsq_residual))
  {
    SLEQP_CALL(sleqp_sparse_vector_clear(lsq_data->lsq_residual));

    SLEQP_CALL(lsq_data->callbacks.lsq_eval(func,
                                            lsq_data->lsq_residual,
                                            lsq_data->func_data));

    lsq_data->has_lsq_residual = true;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE lsq_func_val(SleqpFunc* func,
                                  double* func_val,
                                  void* func_data)
{
  SleqpLSQData* lsq_data = (SleqpLSQData*) func_data;

  *func_val = 0.;

  if(lsq_data->callbacks.additional_func_val)
  {
    //lsq_data->func_data
    //SLEQP_CALL(sleqp_func_val(lsq_data->func, &additional_func_val));
    SLEQP_CALL(lsq_data->callbacks.additional_func_val(func,
                                                       func_val,
                                                       lsq_data->func_data));
  }

  SLEQP_CALL(compute_lsq_residual(func, lsq_data));

  *func_val += .5 * sleqp_sparse_vector_norm_sq(lsq_data->lsq_residual);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
lsq_func_grad(SleqpFunc* func,
              SleqpSparseVec* func_grad,
              void* func_data)
{
  SleqpLSQData* lsq_data = (SleqpLSQData*) func_data;

  SLEQP_CALL(sleqp_sparse_vector_clear(func_grad));

  SLEQP_CALL(compute_lsq_residual(func, lsq_data));

  if(lsq_data->callbacks.additional_func_grad)
  {
    SLEQP_CALL(lsq_data->callbacks.lsq_jac_adjoint(func,
                                                   lsq_data->lsq_residual,
                                                   lsq_data->lsq_grad,
                                                   lsq_data->func_data));

    SLEQP_CALL(lsq_data->callbacks.additional_func_grad(func,
                                                        lsq_data->grad_cache,
                                                        lsq_data->func_data));

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(lsq_data->grad_cache,
                                              lsq_data->lsq_grad,
                                              1.,
                                              1.,
                                              lsq_data->eps,
                                              func_grad));

  }
  else
  {
    SLEQP_CALL(lsq_data->callbacks.lsq_jac_adjoint(func,
                                                   lsq_data->lsq_residual,
                                                   func_grad,
                                                   lsq_data->func_data));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
lsq_func_cons_val(SleqpFunc* func,
                  const SleqpSparseVec* cons_indices,
                  SleqpSparseVec* cons_val,
                  void* func_data)
{
  SleqpLSQData* lsq_data = (SleqpLSQData*) func_data;

  SLEQP_CALL(sleqp_sparse_vector_clear(cons_val));

  if(lsq_data->callbacks.additional_cons_val)
  {
    SLEQP_CALL(lsq_data->callbacks.additional_cons_val(func,
                                                       cons_indices,
                                                       cons_val,
                                                       lsq_data->func_data));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
lsq_func_cons_jac(SleqpFunc* func,
                  const SleqpSparseVec* cons_indices,
                  SleqpSparseMatrix* cons_jac,
                  void* func_data)
{
  SleqpLSQData* lsq_data = (SleqpLSQData*) func_data;

  SLEQP_CALL(sleqp_sparse_matrix_clear(cons_jac));

  if(lsq_data->callbacks.additional_cons_jac)
  {
    SLEQP_CALL(lsq_data->callbacks.additional_cons_jac(func,
                                                       cons_indices,
                                                       cons_jac,
                                                       lsq_data->func_data));
  }

  return SLEQP_OKAY;
}


static SLEQP_RETCODE lsq_func_hess_product(SleqpFunc* func,
                                           const double* func_dual,
                                           const SleqpSparseVec* direction,
                                           const SleqpSparseVec* cons_duals,
                                           SleqpSparseVec* product,
                                           void* func_data)
{
  SleqpLSQData* lsq_data = (SleqpLSQData*) func_data;

  SLEQP_CALL(sleqp_sparse_vector_clear(lsq_data->lsq_forward));
  SLEQP_CALL(sleqp_sparse_vector_clear(lsq_data->lsq_hess_prod));

  const bool additional_term = (lsq_data->levenberg_marquardt != 0.);

  SleqpSparseVec* initial_product_dest = product;

  if(additional_term)
  {
    initial_product_dest = lsq_data->combined_hess_prod;
  }

  if(func_dual)
  {
    if(lsq_data->callbacks.additional_hess_prod)
    {
      SLEQP_CALL(lsq_data->callbacks.additional_hess_prod(func,
                                                          func_dual,
                                                          direction,
                                                          cons_duals,
                                                          lsq_data->hess_prod,
                                                          lsq_data->func_data));

      SLEQP_CALL(lsq_data->callbacks.lsq_jac_forward(func,
                                                     direction,
                                                     lsq_data->lsq_forward,
                                                     lsq_data->func_data));

      SLEQP_CALL(lsq_data->callbacks.lsq_jac_adjoint(func,
                                                     lsq_data->lsq_forward,
                                                     lsq_data->lsq_hess_prod,
                                                     lsq_data->func_data));

      SLEQP_CALL(sleqp_sparse_vector_add_scaled(lsq_data->hess_prod,
                                                lsq_data->lsq_hess_prod,
                                                1.,
                                                (*func_dual),
                                                lsq_data->eps,
                                                initial_product_dest));
    }
    else
    {
      SLEQP_CALL(lsq_data->callbacks.lsq_jac_forward(func,
                                                     direction,
                                                     lsq_data->lsq_forward,
                                                     lsq_data->func_data));

      SLEQP_CALL(lsq_data->callbacks.lsq_jac_adjoint(func,
                                                     lsq_data->lsq_forward,
                                                     initial_product_dest,
                                                     lsq_data->func_data));
    }
  }
  else if(lsq_data->callbacks.additional_hess_prod)
  {
    SLEQP_CALL(lsq_data->callbacks.additional_hess_prod(func,
                                                        func_dual,
                                                        direction,
                                                        cons_duals,
                                                        initial_product_dest,
                                                        lsq_data->func_data));
  }


  if(additional_term)
  {
    SLEQP_CALL(sleqp_sparse_vector_add_scaled(direction,
                                              initial_product_dest,
                                              lsq_data->levenberg_marquardt,
                                              1.,
                                              lsq_data->eps,
                                              product));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE lsq_func_free(void* func_data)
{
  if(!func_data)
  {
    return SLEQP_OKAY;
  }
  SleqpLSQData* lsq_data = (SleqpLSQData*) func_data;

  if(lsq_data->callbacks.func_free)
  {
    SLEQP_CALL(lsq_data->callbacks.func_free(lsq_data->func_data));
  }

  SLEQP_CALL(sleqp_sparse_vector_free(&lsq_data->combined_hess_prod));

  SLEQP_CALL(sleqp_sparse_vector_free(&lsq_data->hess_prod));

  SLEQP_CALL(sleqp_sparse_vector_free(&lsq_data->grad_cache));

  SLEQP_CALL(sleqp_sparse_vector_free(&lsq_data->lsq_grad));

  SLEQP_CALL(sleqp_sparse_vector_free(&lsq_data->lsq_hess_prod));

  SLEQP_CALL(sleqp_sparse_vector_free(&lsq_data->lsq_residual));

  SLEQP_CALL(sleqp_sparse_vector_free(&lsq_data->lsq_forward));

  sleqp_free(&lsq_data);

  return SLEQP_OKAY;
}


SLEQP_RETCODE sleqp_lsq_func_create(SleqpFunc** fstar,
                                    SleqpLSQCallbacks* callbacks,
                                    int num_variables,
                                    int num_constraints,
                                    int num_residuals,
                                    double levenberg_marquardt,
                                    SleqpParams* params,
                                    void* func_data)
{
  SleqpLSQData* data = NULL;

  SLEQP_CALL(sleqp_malloc(&data));

  *data = (SleqpLSQData) {0};

  data->num_variables = num_variables;
  data->num_residuals = num_residuals;

  data->has_lsq_residual = false;

  data->callbacks = *callbacks;

  data->levenberg_marquardt = levenberg_marquardt;

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->lsq_forward,
                                              num_residuals));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->lsq_residual,
                                              num_residuals));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->lsq_hess_prod,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->lsq_grad,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->grad_cache,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->hess_prod,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->combined_hess_prod,
                                              num_variables));

  data->eps = sleqp_params_get(params, SLEQP_PARAM_EPS);

  data->func_data = func_data;

  SleqpFuncCallbacks func_callbacks = {
    .set_value = lsq_func_set_value,
    .func_val  = lsq_func_val,
    .func_grad = lsq_func_grad,
    .cons_val  = lsq_func_cons_val,
    .cons_jac  = lsq_func_cons_jac,
    .hess_prod = lsq_func_hess_product,
    .func_free = lsq_func_free
  };

  SLEQP_CALL(sleqp_func_create(fstar,
                               &func_callbacks,
                               num_variables,
                               num_constraints,
                               data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_lsq_func_set_callbacks(SleqpFunc* func,
                                           SleqpLSQCallbacks* callbacks)
{
  void* func_data = sleqp_func_get_data(func);

  if(!func_data)
  {
    return SLEQP_OKAY;
  }

  SleqpLSQData* lsq_data = (SleqpLSQData*) func_data;

  lsq_data->callbacks = *callbacks;

  return SLEQP_OKAY;
}
