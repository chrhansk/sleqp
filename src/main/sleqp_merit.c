#include "sleqp_merit.h"

#include "sleqp_cmp.h"
#include "sleqp_feas.h"
#include "sleqp_mem.h"
#include "sleqp_util.h"

struct SleqpMeritData
{
  int refcount;

  SleqpProblem* problem;
  SleqpParams* params;

  double* dense_cache;
  int cache_size;

  SleqpSparseVec* jac_dot_sparse;
  SleqpSparseVec* combined_cons_val;

  SleqpSparseVec* multipliers;
  SleqpSparseVec* sparse_cache;

  SleqpFunc* func;
};

SLEQP_RETCODE sleqp_merit_data_create(SleqpMeritData** star,
                                      SleqpProblem* problem,
                                      SleqpParams* params)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpMeritData* merit_data = *star;

  merit_data->refcount = 1;

  merit_data->problem = problem;
  merit_data->func = problem->func;

  SLEQP_CALL(sleqp_params_capture(params));
  merit_data->params = params;

  merit_data->cache_size = SLEQP_MAX(problem->num_constraints,
                                     problem->num_variables);

  SLEQP_CALL(sleqp_calloc(&merit_data->dense_cache, merit_data->cache_size));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&merit_data->jac_dot_sparse,
                                              problem->num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&merit_data->combined_cons_val,
                                              problem->num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&merit_data->multipliers,
                                              problem->num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&merit_data->sparse_cache,
                                              problem->num_constraints));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_merit_func(SleqpMeritData* merit_data,
                               SleqpIterate* iterate,
                               double penalty_parameter,
                               double* merit_value)
{
  *merit_value = sleqp_iterate_get_func_val(iterate);

  SleqpSparseVec* lb = merit_data->problem->cons_lb;
  SleqpSparseVec* ub = merit_data->problem->cons_ub;
  SleqpSparseVec* c = sleqp_iterate_get_cons_val(iterate);

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
                                 const SleqpSparseVec* direction,
                                 double penalty_parameter,
                                 double* merit_value)
{
  SleqpProblem* problem = merit_data->problem;

  SleqpSparseVec* lb = problem->cons_lb;
  SleqpSparseVec* ub = problem->cons_ub;

  const double zero_eps = sleqp_params_get_zero_eps(merit_data->params);

  double objective_dot;

  SLEQP_CALL(sleqp_sparse_vector_dot(sleqp_iterate_get_func_grad(iterate),
                                     direction,
                                     &objective_dot));

  (*merit_value) = sleqp_iterate_get_func_val(iterate) + objective_dot;

  SLEQP_CALL(sleqp_sparse_matrix_vector_product(sleqp_iterate_get_cons_jac(iterate),
                                                direction,
                                                merit_data->dense_cache));

  SLEQP_CALL(sleqp_sparse_vector_from_raw(merit_data->jac_dot_sparse,
                                          merit_data->dense_cache,
                                          problem->num_constraints,
                                          zero_eps));

  SLEQP_CALL(sleqp_sparse_vector_add(merit_data->jac_dot_sparse,
                                     sleqp_iterate_get_cons_val(iterate),
                                     zero_eps,
                                     merit_data->combined_cons_val));

  double total_violation;

  SLEQP_CALL(sleqp_violation_one_norm(problem,
                                      merit_data->combined_cons_val,
                                      zero_eps,
                                      &total_violation));

  (*merit_value) += penalty_parameter * total_violation;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_merit_quadratic(SleqpMeritData* merit_data,
                                    SleqpIterate* iterate,
                                    const double* func_dual,
                                    const SleqpSparseVec* direction,
                                    const SleqpSparseVec* cons_duals,
                                    double penalty_parameter,
                                    double* merit_value)
{
  double linear_merit_value;

  SLEQP_CALL(sleqp_merit_linear(merit_data,
                                iterate,
                                direction,
                                penalty_parameter,
                                &linear_merit_value));

  SleqpFunc* func = merit_data->func;

  double bilinear_product;

  SLEQP_CALL(sleqp_func_hess_bilinear(func,
                                      func_dual,
                                      direction,
                                      cons_duals,
                                      &bilinear_product));

  (*merit_value) = linear_merit_value + (0.5 * bilinear_product);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE merit_data_free(SleqpMeritData** star)
{
  SleqpMeritData* merit_data = *star;

  if(!merit_data)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_sparse_vector_free(&merit_data->sparse_cache));

  SLEQP_CALL(sleqp_sparse_vector_free(&merit_data->multipliers));

  SLEQP_CALL(sleqp_sparse_vector_free(&merit_data->combined_cons_val));

  SLEQP_CALL(sleqp_sparse_vector_free(&merit_data->jac_dot_sparse));

  sleqp_free(&merit_data->dense_cache);

  SLEQP_CALL(sleqp_params_release(&merit_data->params));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_merit_data_capture(SleqpMeritData* merit_data)
{
  ++merit_data->refcount;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_merit_data_release(SleqpMeritData** star)
{
  SleqpMeritData* merit_data = *star;

  if(!merit_data)
  {
    return SLEQP_OKAY;
  }

  if(--merit_data->refcount == 0)
  {
    SLEQP_CALL(merit_data_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
