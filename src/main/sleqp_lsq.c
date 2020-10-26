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

  // num_variables
  SleqpSparseVec* lsq_hess_prod;
  SleqpSparseVec* lsq_grad;
  SleqpSparseVec* grad_cache;
  SleqpSparseVec* hess_prod;

  SleqpSparseVec* combined_hess_prod;

  double eps;

  void* func_data;

} SleqpLSQData;


static SLEQP_RETCODE lsq_func_set_value(SleqpSparseVec* x,
                                        SLEQP_VALUE_REASON reason,
                                        int num_variables,
                                        int* func_grad_nnz,
                                        int* cons_val_nnz,
                                        int* cons_jac_nnz,
                                        void* func_data)
{
  SleqpLSQData* lsq_data = (SleqpLSQData*) func_data;

  SLEQP_CALL(lsq_data->callbacks.set_value(x,
                                           reason,
                                           num_variables,
                                           func_grad_nnz,
                                           cons_val_nnz,
                                           cons_jac_nnz,
                                           lsq_data->func_data));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE lsq_func_eval(int num_variables,
                                   const SleqpSparseVec* cons_indices,
                                   double* func_val,
                                   SleqpSparseVec* func_grad,
                                   SleqpSparseVec* cons_val,
                                   SleqpSparseMatrix* cons_jac,
                                   void* func_data)
{

  SleqpLSQData* lsq_data = (SleqpLSQData*) func_data;

  SLEQP_CALL(sleqp_sparse_vector_clear(lsq_data->lsq_residual));
  SLEQP_CALL(sleqp_sparse_vector_clear(lsq_data->lsq_grad));

  if(func_val)
  {
    *func_val = 0.;
  }

  if(lsq_data->callbacks.eval_additional)
  {
    SleqpSparseVec* grad = NULL;

    if(func_grad)
    {
      grad = lsq_data->grad_cache;
    }

    SLEQP_CALL(lsq_data->callbacks.eval_additional(num_variables,
                                                   cons_indices,
                                                   func_val,
                                                   grad,
                                                   cons_val,
                                                   cons_jac,
                                                   lsq_data->func_data));
  }

  if(func_val || func_grad)
  {
    SLEQP_CALL(sleqp_sparse_vector_clear(lsq_data->lsq_residual));

    SLEQP_CALL(lsq_data->callbacks.lsq_eval(num_variables,
                                            lsq_data->lsq_residual,
                                            lsq_data->func_data));
  }

  if(func_val)
  {
    *func_val += .5 * sleqp_sparse_vector_norm_sq(lsq_data->lsq_residual);
  }

  if(func_grad)
  {
    SLEQP_CALL(lsq_data->callbacks.lsq_jac_adjoint(num_variables,
                                                   lsq_data->lsq_residual,
                                                   lsq_data->lsq_grad,
                                                   lsq_data->func_data));

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(lsq_data->grad_cache,
                                              lsq_data->lsq_grad,
                                              1.,
                                              1.,
                                              lsq_data->eps,
                                              func_grad));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE lsq_func_hess_product(int num_variables,
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
    if(lsq_data->callbacks.hess_prod_additional)
    {
      SLEQP_CALL(lsq_data->callbacks.hess_prod_additional(num_variables,
                                                          func_dual,
                                                          direction,
                                                          cons_duals,
                                                          lsq_data->hess_prod,
                                                          lsq_data->func_data));

      SLEQP_CALL(lsq_data->callbacks.lsq_jac_forward(num_variables,
                                                     direction,
                                                     lsq_data->lsq_forward,
                                                     lsq_data->func_data));

      SLEQP_CALL(lsq_data->callbacks.lsq_jac_adjoint(num_variables,
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
      SLEQP_CALL(lsq_data->callbacks.lsq_jac_forward(num_variables,
                                                     direction,
                                                     lsq_data->lsq_forward,
                                                     lsq_data->func_data));

      SLEQP_CALL(lsq_data->callbacks.lsq_jac_adjoint(num_variables,
                                                     lsq_data->lsq_forward,
                                                     initial_product_dest,
                                                     lsq_data->func_data));
    }
  }
  else if(lsq_data->callbacks.hess_prod_additional)
  {
    SLEQP_CALL(lsq_data->callbacks.hess_prod_additional(num_variables,
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

  sleqp_sparse_vector_free(&lsq_data->lsq_forward);
  sleqp_sparse_vector_free(&lsq_data->lsq_residual);

  sleqp_sparse_vector_free(&lsq_data->lsq_hess_prod);
  sleqp_sparse_vector_free(&lsq_data->lsq_grad);
  sleqp_sparse_vector_free(&lsq_data->grad_cache);

  sleqp_sparse_vector_free(&lsq_data->hess_prod);
  sleqp_sparse_vector_free(&lsq_data->combined_hess_prod);

  sleqp_free(&lsq_data);

  return SLEQP_OKAY;
}


SLEQP_RETCODE sleqp_lsq_func_create(SleqpFunc** fstar,
                                    SleqpLSQCallbacks* callbacks,
                                    int num_variables,
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

  data->callbacks = *callbacks;

  data->levenberg_marquardt = levenberg_marquardt;

  SLEQP_CALL(sleqp_sparse_vector_create(&data->lsq_forward,
                                        num_residuals,
                                        num_residuals));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->lsq_residual,
                                        num_residuals,
                                        num_residuals));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->lsq_hess_prod,
                                        num_variables,
                                        num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->lsq_grad,
                                        num_variables,
                                        num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->lsq_grad,
                                        num_variables,
                                        num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->grad_cache,
                                        num_variables,
                                        num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->hess_prod,
                                        num_variables,
                                        num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->combined_hess_prod,
                                        num_variables,
                                        num_variables));

  data->eps = sleqp_params_get_eps(params);

  data->func_data = func_data;

  SleqpFuncCallbacks func_callbacks = {
    .set_value = lsq_func_set_value,
    .func_eval = lsq_func_eval,
    .hess_prod = lsq_func_hess_product,
    .func_free = lsq_func_free
  };

  SLEQP_CALL(sleqp_func_create(fstar,
                               &func_callbacks,
                               num_variables,
                               data));

  return SLEQP_OKAY;
}
