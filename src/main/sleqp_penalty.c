#include "sleqp_penalty.h"

#include "sleqp_cmp.h"
#include "sleqp_mem.h"

struct SleqpPenalty
{
  SleqpProblem* problem;

  double* jac_dot_dense;
  SleqpSparseVec* jac_dot_sparse;
  SleqpSparseVec* lin_jac_vals;

  SleqpFunc* func;
};

SLEQP_RETCODE sleqp_penalty_create(SleqpPenalty** star,
                                   SleqpProblem* problem,
                                   SleqpFunc* func)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpPenalty* penalty = *star;

  penalty->problem = problem;
  penalty->func = func;

  SLEQP_CALL(sleqp_calloc(&penalty->jac_dot_dense, problem->num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create(&penalty->jac_dot_sparse,
                                        problem->num_constraints,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&penalty->lin_jac_vals,
                                        problem->num_constraints,
                                        0));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_penalty_func(SleqpPenalty* penalty,
                                 SleqpIterate* iterate,
                                 double penalty_parameter,
                                 double* penalty_value)
{
  *penalty_value = iterate->func_val;

  SleqpSparseVec* lb = penalty->problem->cons_lb;
  SleqpSparseVec* ub = penalty->problem->cons_ub;
  SleqpSparseVec* c = iterate->cons_val;

  size_t k_c = 0, k_lb = 0, k_ub = 0;

  while(k_lb < lb->nnz || k_ub < ub->nnz || k_c < c->nnz)
  {
    size_t i = c->dim + 1;

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

  SLEQP_CALL(sleqp_sparse_matrix_vector_product(iterate->cons_jac,
                                                direction,
                                                penalty_data->jac_dot_dense));

  SLEQP_CALL(sleqp_sparse_vector_from_raw(penalty_data->jac_dot_sparse,
                                          penalty_data->jac_dot_dense,
                                          problem->num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_add(penalty_data->jac_dot_sparse,
                                     iterate->cons_val,
                                     1., 1.,
                                     lin));

  size_t k_l = 0, k_lb = 0, k_ub = 0;

  while(k_l < lin->nnz || k_lb < lb->nnz || k_ub < ub->nnz)
  {
    size_t i = lin->dim + 1;

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

SLEQP_RETCODE sleqp_penalty_quadratic(SleqpPenalty* penalty_data,
                                      SleqpIterate* iterate,
                                      SleqpSparseVec* direction,
                                      double penalty_parameter,
                                      double* penalty_value)
{
  double linear_penalty_value;
  double func_val = 1.;

  SleqpFunc* func = penalty_data->problem->func;

  SLEQP_CALL(sleqp_penalty_linear(penalty_data,
                                  iterate,
                                  direction,
                                  penalty_parameter,
                                  &linear_penalty_value));

  SLEQP_CALL(sleqp_hess_eval_bilinear(func,
                                      &func_val,
                                      direction,
                                      iterate->cons_dual,
                                      penalty_value));

  *penalty_value *= 0.5;
  *penalty_value += linear_penalty_value;


  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_penalty_free(SleqpPenalty** star)
{
  SleqpPenalty* penalty = *star;

  SLEQP_CALL(sleqp_sparse_vector_free(&penalty->lin_jac_vals));

  SLEQP_CALL(sleqp_sparse_vector_free(&penalty->jac_dot_sparse));

  sleqp_free(&penalty->jac_dot_dense);

  sleqp_free(star);

  return SLEQP_OKAY;
}
