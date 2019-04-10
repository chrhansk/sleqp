#include "sleqp_scale.h"

#include <math.h>

#include "sleqp_mem.h"

struct SleqpScalingData
{
  SleqpProblem* problem;

  int func_scale;
  int* var_scales;
  int* cons_scales;

  SleqpSparseVec* unscaled_value;

  SleqpSparseVec* scaled_direction;
  SleqpSparseVec* scaled_cons_duals;

  SleqpFunc* scaled_func;
  SleqpProblem* scaled_problem;
};

static SLEQP_RETCODE apply_const_scaling(SleqpSparseVec* vec,
                                         int scale)
{
  if(scale == 0)
  {
    return SLEQP_OKAY;
  }

  for(int k = 0; k < vec->nnz; ++k)
  {
    vec->data[k] = ldexp(vec->data[k], scale);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE apply_scaling(SleqpSparseVec* vec,
                                   int* scales,
                                   int offset)
{
  for(int k = 0; k < vec->nnz; ++k)
  {
    vec->data[k] = ldexp(vec->data[k], scales[vec->indices[k]] + offset);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE apply_unscaling(SleqpSparseVec* vec,
                                     int* scales,
                                     int offset)
{
  for(int k = 0; k < vec->nnz; ++k)
  {
    vec->data[k] = ldexp(vec->data[k], -1*scales[vec->indices[k]] + offset);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE apply_quad_unscaling(SleqpSparseVec* vec,
                                          int* scales)
{
  for(int k = 0; k < vec->nnz; ++k)
  {
    vec->data[k] = ldexp(vec->data[k], -2*scales[vec->indices[k]]);
  }

  return SLEQP_OKAY;
}


static SLEQP_RETCODE
scaled_func_set_value(SleqpSparseVec* scaled_value,
                      int num_variables,
                      int* func_grad_nnz,
                      int* cons_val_nnz,
                      int* cons_jac_nnz,
                      void* func_data)
{
  SleqpScalingData* scaling_data = (SleqpScalingData*) func_data;

  SLEQP_CALL(sleqp_sparse_vector_copy(scaled_value,
                                      scaling_data->unscaled_value));

  SLEQP_CALL(apply_unscaling(scaling_data->unscaled_value,
                             scaling_data->var_scales,
                             0));

  SLEQP_CALL(sleqp_func_set_value(scaling_data->problem->func,
                                  scaling_data->unscaled_value,
                                  func_grad_nnz,
                                  cons_val_nnz,
                                  cons_jac_nnz));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
scaled_func_eval(int num_variables,
                 SleqpSparseVec* cons_indices,
                 double* func_val,
                 SleqpSparseVec* func_grad,
                 SleqpSparseVec* cons_val,
                 SleqpSparseMatrix* cons_jac,
                 void* func_data)
{
  SleqpScalingData* scaling_data = (SleqpScalingData*) func_data;

  SLEQP_CALL(sleqp_func_eval(scaling_data->problem->func,
                             cons_indices,
                             func_val,
                             func_grad,
                             cons_val,
                             cons_jac));

  if(func_val)
  {
    (*func_val) = sleqp_scale_func_val(scaling_data, (*func_val));
  }

  if(func_grad)
  {
    SLEQP_CALL(sleqp_scale_func_grad(scaling_data, func_grad));
  }

  if(cons_val)
  {
    SLEQP_CALL(sleqp_scale_cons_val(scaling_data, cons_val));
  }

  if(cons_jac)
  {
    SLEQP_CALL(sleqp_scale_cons_jac(scaling_data, cons_jac));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
scaled_func_hess_prod(int num_variables,
                      double* func_dual,
                      SleqpSparseVec* direction,
                      SleqpSparseVec* cons_duals,
                      SleqpSparseVec* product,
                      void* func_data)
{
  SleqpScalingData* scaling_data = (SleqpScalingData*) func_data;

  SLEQP_CALL(sleqp_sparse_vector_copy(direction,
                                      scaling_data->scaled_direction));

  SLEQP_CALL(apply_quad_unscaling(scaling_data->scaled_direction,
                                  scaling_data->var_scales));

  SLEQP_CALL(sleqp_sparse_vector_copy(cons_duals,
                                      scaling_data->scaled_cons_duals));

  SLEQP_CALL(apply_scaling(scaling_data->scaled_cons_duals,
                           scaling_data->cons_scales,
                           -1. * scaling_data->func_scale));

  SLEQP_CALL(sleqp_func_hess_prod(scaling_data->problem->func,
                                  func_dual,
                                  scaling_data->scaled_direction,
                                  scaling_data->scaled_cons_duals,
                                  product));

  SLEQP_CALL(apply_const_scaling(product, scaling_data->func_scale));


  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_scaling_create(SleqpScalingData** star,
                                   SleqpProblem* problem,
                                   SleqpParams* params)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpScalingData* scaling_data = *star;

  *scaling_data = (SleqpScalingData) {0};

  scaling_data->problem = problem;

  scaling_data->func_scale = 0;

  SLEQP_CALL(sleqp_calloc(&(scaling_data->var_scales),
                          problem->num_variables));

  for(int j = 0; j < problem->num_variables; ++j)
  {
    scaling_data->var_scales[j] = 0;
  }

  SLEQP_CALL(sleqp_calloc(&(scaling_data->cons_scales),
                          problem->num_constraints));

  for(int i = 0; i < problem->num_constraints; ++i)
  {
    scaling_data->cons_scales[i] = 0;
  }

  SLEQP_CALL(sleqp_sparse_vector_create(&(scaling_data->unscaled_value),
                                        problem->num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&(scaling_data->scaled_direction),
                                        problem->num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&(scaling_data->scaled_cons_duals),
                                        problem->num_constraints,
                                        0));

  SLEQP_CALL(sleqp_func_create(&(scaling_data->scaled_func),
                               scaled_func_set_value,
                               scaled_func_eval,
                               scaled_func_hess_prod,
                               problem->num_variables,
                               scaling_data));

  SLEQP_CALL(sleqp_problem_create(&(scaling_data->scaled_problem),
                                  scaling_data->scaled_func,
                                  params,
                                  problem->var_lb,
                                  problem->var_ub,
                                  problem->cons_lb,
                                  problem->cons_ub));

  return SLEQP_OKAY;
}

SleqpProblem* sleqp_scaling_get_scaled_problem(SleqpScalingData* scaling)
{
  return scaling->scaled_problem;
}

SLEQP_RETCODE sleqp_scaling_set_func_scale(SleqpScalingData* scaling,
                                           int scale)
{
  scaling->func_scale = scale;
  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_scaling_set_var_scale(SleqpScalingData* scaling,
                                          int index,
                                          int scale)
{
  if(index < 0 || index >= scaling->problem->num_variables)
  {
    return SLEQP_ILLEGAL_ARGUMENT;
  }

  scaling->cons_scales[index] = scale;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_scaling_set_cons_scale(SleqpScalingData* scaling,
                                           int index,
                                           int scale)
{
  if(index < 0 || index >= scaling->problem->num_constraints)
  {
    return SLEQP_ILLEGAL_ARGUMENT;
  }

  scaling->var_scales[index] = scale;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_scaling_flush(SleqpScalingData* scaling)
{
  SleqpProblem* problem = scaling->problem;
  SleqpProblem* scaled_problem = scaling->scaled_problem;

  SLEQP_CALL(sleqp_sparse_vector_copy(problem->var_lb,
                                      scaled_problem->var_lb));

  SLEQP_CALL(apply_scaling(scaled_problem->var_lb,
                           scaling->var_scales,
                           0));

  SLEQP_CALL(sleqp_sparse_vector_copy(problem->var_ub,
                                      scaled_problem->var_ub));

  SLEQP_CALL(apply_scaling(scaled_problem->var_ub,
                           scaling->var_scales,
                           0));

  SLEQP_CALL(sleqp_sparse_vector_copy(problem->cons_lb,
                                      scaled_problem->cons_lb));

  SLEQP_CALL(apply_scaling(scaled_problem->cons_lb,
                           scaling->cons_scales,
                           0));

  SLEQP_CALL(sleqp_sparse_vector_copy(problem->cons_ub,
                                      scaled_problem->cons_ub));

  SLEQP_CALL(apply_scaling(scaled_problem->cons_ub,
                           scaling->cons_scales,
                           0));

  return SLEQP_OKAY;
}

double sleqp_unscale_func_val(SleqpScalingData* scaling,
                              double scaled_func_val)
{
  return ldexp(scaled_func_val, -1*scaling->func_scale);
}

double sleqp_scale_func_val(SleqpScalingData* scaling,
                            double func_val)
{
  return ldexp(func_val, scaling->func_scale);
}

SLEQP_RETCODE sleqp_scale_point(SleqpScalingData* scaling,
                                  SleqpSparseVec* point)
{
  SLEQP_CALL(apply_scaling(point, scaling->var_scales, 0));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_scale_func_grad(SleqpScalingData* scaling,
                                    SleqpSparseVec* func_grad)
{
  SLEQP_CALL(apply_unscaling(func_grad,
                             scaling->var_scales,
                             scaling->func_scale));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_scale_cons_val(SleqpScalingData* scaling,
                                   SleqpSparseVec* cons_val)
{
  SLEQP_CALL(apply_scaling(cons_val, scaling->cons_scales, 0));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_scale_cons_jac(SleqpScalingData* scaling,
                                   SleqpSparseMatrix* cons_jac)
{
  int col = 0;

  for(int index = 0; index < cons_jac->nnz; ++index)
  {
    while(index >= cons_jac->cols[col + 1])
    {
      ++col;
    }

    const int row = cons_jac->rows[index];

    cons_jac->data[index] = ldexp(cons_jac->data[index],
                                  scaling->cons_scales[row] -
                                  scaling->var_scales[col]);
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_scale_cons_duals(SleqpScalingData* scaling,
                                     SleqpSparseVec* cons_duals)
{
  SLEQP_CALL(apply_unscaling(cons_duals,
                             scaling->cons_scales,
                             -1* scaling->func_scale));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_scale_var_duals(SleqpScalingData* scaling,
                                      SleqpSparseVec* var_duals)
{
  SLEQP_CALL(apply_unscaling(var_duals,
                             scaling->var_scales,
                             -1* scaling->func_scale));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_scale_iterate(SleqpScalingData* scaling,
                                    SleqpIterate* scaled_iterate)
{
  SLEQP_CALL(sleqp_scale_point(scaling, scaled_iterate->x));

  scaled_iterate->func_val = sleqp_scale_func_val(scaling, scaled_iterate->func_val);

  SLEQP_CALL(sleqp_scale_func_grad(scaling, scaled_iterate->func_grad));

  SLEQP_CALL(sleqp_scale_cons_val(scaling, scaled_iterate->cons_val));

  SLEQP_CALL(sleqp_scale_cons_jac(scaling, scaled_iterate->cons_jac));

  SLEQP_CALL(sleqp_scale_cons_duals(scaling, scaled_iterate->cons_dual));

  SLEQP_CALL(sleqp_scale_var_duals(scaling, scaled_iterate->vars_dual));

  return SLEQP_OKAY;
}




SLEQP_RETCODE sleqp_unscale_point(SleqpScalingData* scaling,
                                  SleqpSparseVec* scaled_point)
{
  SLEQP_CALL(apply_unscaling(scaled_point, scaling->var_scales, 0));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_unscale_func_grad(SleqpScalingData* scaling,
                                      SleqpSparseVec* scaled_func_grad)
{
  SLEQP_CALL(apply_scaling(scaled_func_grad,
                           scaling->var_scales,
                           (-1)*scaling->func_scale));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_unscale_cons_val(SleqpScalingData* scaling,
                                     SleqpSparseVec* scaled_cons_val)
{
  SLEQP_CALL(apply_unscaling(scaled_cons_val, scaling->cons_scales, 0));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_unscale_cons_jac(SleqpScalingData* scaling,
                                     SleqpSparseMatrix* scaled_cons_jac)
{
  int col = 0;

  for(int index = 0; index < scaled_cons_jac->nnz; ++index)
  {
    while(index >= scaled_cons_jac->cols[col + 1])
    {
      ++col;
    }

    const int row = scaled_cons_jac->rows[index];

    scaled_cons_jac->data[index] = ldexp(scaled_cons_jac->data[index],
                                         scaling->var_scales[col] -
                                         scaling->cons_scales[row]);
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_unscale_cons_duals(SleqpScalingData* scaling,
                                       SleqpSparseVec* scaled_cons_duals)
{
  SLEQP_CALL(apply_scaling(scaled_cons_duals,
                           scaling->cons_scales,
                           -1* scaling->func_scale));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_unscale_var_duals(SleqpScalingData* scaling,
                                      SleqpSparseVec* scaled_var_duals)
{
  SLEQP_CALL(apply_scaling(scaled_var_duals,
                           scaling->var_scales,
                           -1* scaling->func_scale));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_unscale_iterate(SleqpScalingData* scaling,
                                    SleqpIterate* scaled_iterate)
{
  SLEQP_CALL(sleqp_unscale_point(scaling, scaled_iterate->x));

  scaled_iterate->func_val = sleqp_unscale_func_val(scaling, scaled_iterate->func_val);

  SLEQP_CALL(sleqp_unscale_func_grad(scaling, scaled_iterate->func_grad));

  SLEQP_CALL(sleqp_unscale_cons_val(scaling, scaled_iterate->cons_val));

  SLEQP_CALL(sleqp_unscale_cons_jac(scaling, scaled_iterate->cons_jac));

  SLEQP_CALL(sleqp_unscale_cons_duals(scaling, scaled_iterate->cons_dual));

  SLEQP_CALL(sleqp_unscale_var_duals(scaling, scaled_iterate->vars_dual));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_scaling_free(SleqpScalingData** star)
{
  SleqpScalingData* scaling_data = *star;

  if(!scaling_data)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_problem_free(&(scaling_data->scaled_problem)));

  SLEQP_CALL(sleqp_func_free(&(scaling_data->scaled_func)));

  SLEQP_CALL(sleqp_sparse_vector_free(&(scaling_data->scaled_cons_duals)));

  SLEQP_CALL(sleqp_sparse_vector_free(&(scaling_data->scaled_direction)));

  SLEQP_CALL(sleqp_sparse_vector_free(&(scaling_data->unscaled_value)));

  sleqp_free(&(scaling_data->cons_scales));

  sleqp_free(&(scaling_data->var_scales));

  sleqp_free(star);

  *star = NULL;

  return SLEQP_OKAY;
}
