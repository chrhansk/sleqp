#include "sleqp_merit.h"

#include "sleqp_cmp.h"
#include "sleqp_mem.h"

struct SleqpMeritData
{
  SleqpProblem* problem;
  SleqpParams* params;

  double* dense_cache;
  int cache_size;

  SleqpSparseVec* jac_dot_sparse;
  SleqpSparseVec* lin_jac_vals;

  SleqpSparseVec* multipliers;
  SleqpSparseVec* sparse_cache;

  SleqpFunc* func;
};

static void reset_cache(double* values, int size)
{
  for(int i = 0; i < size; ++i)
  {
    values[i] = 0.;
  }
}

SLEQP_RETCODE sleqp_merit_data_create(SleqpMeritData** star,
                                      SleqpProblem* problem,
                                      SleqpParams* params)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpMeritData* merit_data = *star;

  merit_data->problem = problem;
  merit_data->func = problem->func;
  merit_data->params = params;

  merit_data->cache_size = SLEQP_MAX(problem->num_constraints,
                                     problem->num_variables);

  SLEQP_CALL(sleqp_calloc(&merit_data->dense_cache, merit_data->cache_size));

  SLEQP_CALL(sleqp_sparse_vector_create(&merit_data->jac_dot_sparse,
                                        problem->num_constraints,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&merit_data->lin_jac_vals,
                                        problem->num_constraints,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&merit_data->multipliers,
                                        problem->num_constraints,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&merit_data->sparse_cache,
                                        problem->num_constraints,
                                        0));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_merit_func(SleqpMeritData* merit_data,
                               SleqpIterate* iterate,
                               double penalty_parameter,
                               double* merit_value)
{
  *merit_value = iterate->func_val;

  SleqpSparseVec* lb = merit_data->problem->cons_lb;
  SleqpSparseVec* ub = merit_data->problem->cons_ub;
  SleqpSparseVec* c = iterate->cons_val;

  int k_c = 0, k_lb = 0, k_ub = 0;

  while(k_lb < lb->nnz || k_ub < ub->nnz || k_c < c->nnz)
  {
    int i = c->dim + 1;

    bool valid_lb = (k_lb < lb->nnz);
    bool valid_ub = (k_ub < ub->nnz);
    bool valid_c = (k_c < c->nnz);

    if(valid_lb)
    {
      i = SLEQP_MIN(i, lb->indices[k_lb]);
    }

    if(valid_ub)
    {
      i = SLEQP_MIN(i, ub->indices[k_ub]);
    }

    if(valid_c)
    {
      i = SLEQP_MIN(i, c->indices[k_c]);
    }

    valid_lb = valid_lb && lb->indices[k_lb] == i;
    valid_ub = valid_ub && ub->indices[k_ub] == i;
    valid_c = valid_c && c->indices[k_c] == i;

    double val_lb = valid_lb ? lb->data[k_lb] : 0;
    double val_ub = valid_ub ? ub->data[k_ub] : 0;
    double val_c = valid_c ? c->data[k_c] : 0;

    *merit_value += penalty_parameter * SLEQP_MAX(val_c - val_ub, 0);
    *merit_value += penalty_parameter * SLEQP_MAX(val_lb - val_c, 0);

    if(valid_lb)
    {
      ++k_lb;
    }

    if(valid_ub)
    {
      ++k_ub;
    }

    if(valid_c)
    {
      ++k_c;
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_merit_linear(SleqpMeritData* merit_data,
                                 SleqpIterate* iterate,
                                 SleqpSparseVec* direction,
                                 double penalty_parameter,
                                 double* merit_value)
{
  SleqpProblem* problem = merit_data->problem;

  SleqpSparseVec* lin = merit_data->lin_jac_vals;
  SleqpSparseVec* lb = problem->cons_lb;
  SleqpSparseVec* ub = problem->cons_ub;

  const double zero_eps = sleqp_params_get_zero_eps(merit_data->params);

  SLEQP_CALL(sleqp_sparse_vector_dot(iterate->func_grad,
                                     direction,
                                     merit_value));

  *merit_value += iterate->func_val;

  reset_cache(merit_data->dense_cache, merit_data->cache_size);

  SLEQP_CALL(sleqp_sparse_matrix_vector_product(iterate->cons_jac,
                                                direction,
                                                merit_data->dense_cache));

  SLEQP_CALL(sleqp_sparse_vector_from_raw(merit_data->jac_dot_sparse,
                                          merit_data->dense_cache,
                                          problem->num_constraints,
                                          zero_eps));

  SLEQP_CALL(sleqp_sparse_vector_add(merit_data->jac_dot_sparse,
                                     iterate->cons_val,
                                     zero_eps,
                                     lin));

  int k_l = 0, k_lb = 0, k_ub = 0;

  while(k_l < lin->nnz || k_lb < lb->nnz || k_ub < ub->nnz)
  {
    int i = lin->dim + 1;

    bool valid_lb = (k_lb < lb->nnz);
    bool valid_ub = (k_ub < ub->nnz);
    bool valid_lin = (k_l < lin->nnz);

    if(valid_lb)
    {
      i = SLEQP_MIN(i, lb->indices[k_lb]);
    }

    if(valid_ub)
    {
      i = SLEQP_MIN(i, ub->indices[k_ub]);
    }

    if(valid_lin)
    {
      i = SLEQP_MIN(i, lin->indices[k_l]);
    }

    valid_lb = valid_lb && lb->indices[k_lb] == i;
    valid_ub = valid_ub && ub->indices[k_ub] == i;
    valid_lin = valid_lin && lin->indices[k_l] == i;

    double val_lb = valid_lb ? lb->data[k_lb] : 0;
    double val_ub = valid_ub ? ub->data[k_ub] : 0;
    double val_lin = valid_lin ? lin->data[k_l] : 0;


    *merit_value += penalty_parameter * SLEQP_MAX(val_lin - val_ub, 0);
    *merit_value += penalty_parameter * SLEQP_MAX(val_lb - val_lin, 0);

    if(valid_lb)
    {
      ++k_lb;
    }

    if(valid_ub)
    {
      ++k_ub;
    }

    if(valid_lin)
    {
      ++k_l;
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_merit_linear_gradient(SleqpMeritData* merit_data,
                                          SleqpIterate* iterate,
                                          SleqpSparseVec* direction,
                                          double penalty_parameter,
                                          SleqpSparseVec* gradient)
{
  SleqpProblem* problem = merit_data->problem;

  const double zero_eps = sleqp_params_get_zero_eps(merit_data->params);

  SleqpSparseVec* x = iterate->primal;
  SleqpSparseVec* cons_vals = iterate->cons_val;

  SLEQP_CALL(sleqp_get_violated_multipliers(problem,
                                            x,
                                            cons_vals,
                                            penalty_parameter,
                                            merit_data->multipliers,
                                            NULL,
                                            zero_eps));

  SLEQP_CALL(sleqp_sparse_vector_clear(merit_data->sparse_cache));

  SLEQP_CALL(sleqp_sparse_vector_resize(merit_data->sparse_cache,
                                        problem->num_variables));

  SLEQP_CALL(sleqp_sparse_matrix_trans_vector_product(iterate->cons_jac,
                                                      merit_data->multipliers,
                                                      zero_eps,
                                                      merit_data->sparse_cache));

  SLEQP_CALL(sleqp_sparse_vector_add(merit_data->sparse_cache,
                                     iterate->func_grad,
                                     zero_eps,
                                     gradient));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_merit_data_free(SleqpMeritData** star)
{
  SleqpMeritData* merit_data = *star;

  if(!merit_data)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_sparse_vector_free(&merit_data->sparse_cache));

  SLEQP_CALL(sleqp_sparse_vector_free(&merit_data->multipliers));

  SLEQP_CALL(sleqp_sparse_vector_free(&merit_data->lin_jac_vals));

  SLEQP_CALL(sleqp_sparse_vector_free(&merit_data->jac_dot_sparse));

  sleqp_free(&merit_data->dense_cache);

  sleqp_free(star);

  return SLEQP_OKAY;
}
