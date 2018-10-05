#include "sleqp_penalty.h"

#include "sleqp_cmp.h"
#include "sleqp_mem.h"

struct SleqpPenalty
{
  SleqpProblem* problem;

  double* dense_cache;
  int cache_size;

  SleqpSparseVec* jac_dot_sparse;
  SleqpSparseVec* lin_jac_vals;

  SleqpFunc* func;
};

static void reset_cache(double* values, int size)
{
  for(int i = 0; i < size; ++i)
  {
    values[i] = 0.;
  }
}

SLEQP_RETCODE sleqp_penalty_create(SleqpPenalty** star,
                                   SleqpProblem* problem,
                                   SleqpFunc* func)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpPenalty* penalty_data = *star;

  penalty_data->problem = problem;
  penalty_data->func = func;

  penalty_data->cache_size = SLEQP_MAX(problem->num_constraints,
                                       problem->num_variables);

  SLEQP_CALL(sleqp_calloc(&penalty_data->dense_cache, penalty_data->cache_size));

  SLEQP_CALL(sleqp_sparse_vector_create(&penalty_data->jac_dot_sparse,
                                        problem->num_constraints,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&penalty_data->lin_jac_vals,
                                        problem->num_constraints,
                                        0));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_penalty_func(SleqpPenalty* penalty_data,
                                 SleqpIterate* iterate,
                                 double penalty_parameter,
                                 double* penalty_value)
{
  *penalty_value = iterate->func_val;

  SleqpSparseVec* lb = penalty_data->problem->cons_lb;
  SleqpSparseVec* ub = penalty_data->problem->cons_ub;
  SleqpSparseVec* c = iterate->cons_val;

  int k_c = 0, k_lb = 0, k_ub = 0;

  while(k_lb < lb->nnz || k_ub < ub->nnz || k_c < c->nnz)
  {
    int i = c->dim + 1;

    SLEQP_Bool valid_lb = (k_lb < lb->nnz);
    SLEQP_Bool valid_ub = (k_ub < ub->nnz);
    SLEQP_Bool valid_c = (k_c < c->nnz);

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

    double val_lb = valid_lb ? lb->data[i] : 0;
    double val_ub = valid_ub ? ub->data[i] : 0;
    double val_c = valid_c ? c->data[i] : 0;

    *penalty_value += penalty_parameter * SLEQP_MAX(val_c - val_ub, 0);
    *penalty_value += penalty_parameter * SLEQP_MAX(val_lb - val_c, 0);

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

SLEQP_RETCODE sleqp_penalty_linear(SleqpPenalty* penalty_data,
                                   SleqpIterate* iterate,
                                   SleqpSparseVec* direction,
                                   double penalty_parameter,
                                   double* penalty_value)
{
  SleqpProblem* problem = penalty_data->problem;

  SleqpSparseVec* lin = penalty_data->lin_jac_vals;
  SleqpSparseVec* lb = problem->cons_lb;
  SleqpSparseVec* ub = problem->cons_ub;

  SLEQP_CALL(sleqp_sparse_vector_dot(iterate->func_grad,
                                     direction,
                                     penalty_value));

  reset_cache(penalty_data->dense_cache, penalty_data->cache_size);

  SLEQP_CALL(sleqp_sparse_matrix_vector_product(iterate->cons_jac,
                                                direction,
                                                penalty_data->dense_cache));

  SLEQP_CALL(sleqp_sparse_vector_from_raw(penalty_data->jac_dot_sparse,
                                          penalty_data->dense_cache,
                                          problem->num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_add(penalty_data->jac_dot_sparse,
                                     iterate->cons_val,
                                     1., 1.,
                                     lin));

  int k_l = 0, k_lb = 0, k_ub = 0;

  while(k_l < lin->nnz || k_lb < lb->nnz || k_ub < ub->nnz)
  {
    int i = lin->dim + 1;

    SLEQP_Bool valid_lb = (k_lb < lb->nnz);
    SLEQP_Bool valid_ub = (k_ub < ub->nnz);
    SLEQP_Bool valid_lin = (k_l < lin->nnz);

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

    double val_lb = valid_lb ? lb->data[i] : 0;
    double val_ub = valid_ub ? ub->data[i] : 0;
    double val_lin = valid_lin ? lin->data[i] : 0;


    *penalty_value += penalty_parameter * SLEQP_MAX(val_lin - val_ub, 0);
    *penalty_value += penalty_parameter * SLEQP_MAX(val_lb - val_lin, 0);

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

SLEQP_RETCODE sleqp_penalty_quadratic_gradient(SleqpPenalty* penalty_data,
                                               SleqpIterate* iterate,
                                               double penalty_parameter,
                                               SleqpSparseVec* gradient)
{
  reset_cache(penalty_data->dense_cache, penalty_data->cache_size);

  SleqpSparseVec* func_grad = iterate->func_grad;
  double* cache = penalty_data->dense_cache;
  SleqpProblem* problem = penalty_data->problem;

  SLEQP_ACTIVE_STATE* cons_states = sleqp_active_set_cons_states(iterate->active_set);
  SleqpSparseMatrix* cons_jac = iterate->cons_jac;

  for(int k = 0; k < func_grad->nnz; ++k)
  {
    cache[func_grad->indices[k]] = func_grad->data[k];
  }

  int col = 0;

  for(int index = 0; index < cons_jac->nnz; ++index)
  {
    while(index >= cons_jac->cols[col + 1])
    {
      ++col;
    }

    switch (cons_states[col]) {
    case SLEQP_ACTIVE_UPPER:
      cache[cons_jac->rows[index]] += penalty_parameter * cons_jac->data[index];
      break;
    case SLEQP_ACTIVE_LOWER:
      cache[cons_jac->rows[index]] += -1. * penalty_parameter * cons_jac->data[index];
      break;
    case SLEQP_INACTIVE:
      break;
    }
  }

  SLEQP_CALL(sleqp_sparse_vector_from_raw(gradient, cache, problem->num_variables));

  return SLEQP_OKAY;
}


SLEQP_RETCODE sleqp_penalty_free(SleqpPenalty** star)
{
  SleqpPenalty* penalty_data = *star;

  SLEQP_CALL(sleqp_sparse_vector_free(&penalty_data->lin_jac_vals));

  SLEQP_CALL(sleqp_sparse_vector_free(&penalty_data->jac_dot_sparse));

  sleqp_free(&penalty_data->dense_cache);

  sleqp_free(star);

  return SLEQP_OKAY;
}
