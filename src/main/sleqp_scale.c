#include "sleqp_scale.h"

#include <math.h>

#include "sleqp_cmp.h"
#include "sleqp_mem.h"

int max_weight = 10;
int min_weight = -10;

#define CLIP_WEIGHT(w) \
  (w) = SLEQP_MIN(SLEQP_MAX((w), min_weight), max_weight)

struct SleqpScalingData
{
  SleqpProblem* problem;
  SleqpParams* params;

  int func_scale;
  int* var_weights;
  int* cons_weights;

  double* min_cache;
  double* max_cache;

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

static SLEQP_RETCODE apply_quad_scaling(SleqpSparseVec* vec,
                                        int* scales)
{
  for(int k = 0; k < vec->nnz; ++k)
  {
    vec->data[k] = ldexp(vec->data[k], 2*scales[vec->indices[k]]);
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
  SleqpScalingData* scaling = (SleqpScalingData*) func_data;

  SLEQP_CALL(sleqp_sparse_vector_copy(scaled_value,
                                      scaling->unscaled_value));

  SLEQP_CALL(sleqp_unscale_point(scaling,
                                 scaling->unscaled_value));

  SLEQP_CALL(sleqp_func_set_value(scaling->problem->func,
                                  scaling->unscaled_value,
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
  SleqpScalingData* scaling = (SleqpScalingData*) func_data;

  SLEQP_CALL(sleqp_func_eval(scaling->problem->func,
                             cons_indices,
                             func_val,
                             func_grad,
                             cons_val,
                             cons_jac));

  if(func_val)
  {
    (*func_val) = sleqp_scale_func_val(scaling, (*func_val));
  }

  if(func_grad)
  {
    SLEQP_CALL(sleqp_scale_func_grad(scaling, func_grad));
  }

  if(cons_val)
  {
    SLEQP_CALL(sleqp_scale_cons_val(scaling, cons_val));
  }

  if(cons_jac)
  {
    SLEQP_CALL(sleqp_scale_cons_jac(scaling, cons_jac));
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
  SleqpScalingData* scaling = (SleqpScalingData*) func_data;

  SLEQP_CALL(sleqp_sparse_vector_copy(direction,
                                      scaling->scaled_direction));

  SLEQP_CALL(apply_quad_scaling(scaling->scaled_direction,
                                scaling->var_weights));

  SLEQP_CALL(sleqp_sparse_vector_copy(cons_duals,
                                      scaling->scaled_cons_duals));

  SLEQP_CALL(apply_scaling(scaling->scaled_cons_duals,
                           scaling->cons_weights,
                           -1. * scaling->func_scale));

  SLEQP_CALL(sleqp_func_hess_prod(scaling->problem->func,
                                  func_dual,
                                  scaling->scaled_direction,
                                  scaling->scaled_cons_duals,
                                  product));

  SLEQP_CALL(apply_const_scaling(product, scaling->func_scale));


  return SLEQP_OKAY;
}

static SLEQP_RETCODE reset_scaling(SleqpScalingData* scaling)
{
  SleqpProblem* problem = scaling->problem;

  for(int j = 0; j < problem->num_variables; ++j)
  {
    scaling->var_weights[j] = 0;
  }

  for(int i = 0; i < problem->num_constraints; ++i)
  {
    scaling->cons_weights[i] = 0;
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_scaling_create(SleqpScalingData** star,
                                   SleqpProblem* problem,
                                   SleqpParams* params)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpScalingData* scaling = *star;

  *scaling = (SleqpScalingData) {0};

  scaling->problem = problem;

  scaling->params = params;

  SLEQP_CALL(sleqp_calloc(&(scaling->var_weights),
                          problem->num_variables));

  SLEQP_CALL(sleqp_calloc(&(scaling->cons_weights),
                          problem->num_constraints));

  scaling->func_scale = 0;

  SLEQP_CALL(reset_scaling(scaling));

  SLEQP_CALL(sleqp_sparse_vector_create(&(scaling->unscaled_value),
                                        problem->num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&(scaling->scaled_direction),
                                        problem->num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&(scaling->scaled_cons_duals),
                                        problem->num_constraints,
                                        0));

  SLEQP_CALL(sleqp_func_create(&(scaling->scaled_func),
                               scaled_func_set_value,
                               scaled_func_eval,
                               scaled_func_hess_prod,
                               problem->num_variables,
                               scaling));

  SLEQP_CALL(sleqp_problem_create(&(scaling->scaled_problem),
                                  scaling->scaled_func,
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

SLEQP_RETCODE sleqp_scaling_set_func_weight(SleqpScalingData* scaling,
                                            int weight)
{
  scaling->func_scale = weight;
  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_scaling_set_var_weight(SleqpScalingData* scaling,
                                           int index,
                                           int weight)
{
  if(index < 0 || index >= scaling->problem->num_variables)
  {
    return SLEQP_ILLEGAL_ARGUMENT;
  }

  scaling->cons_weights[index] = weight;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_scaling_set_cons_weight(SleqpScalingData* scaling,
                                            int index,
                                            int weight)
{
  if(index < 0 || index >= scaling->problem->num_constraints)
  {
    return SLEQP_ILLEGAL_ARGUMENT;
  }

  scaling->var_weights[index] = weight;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_scaling_flush(SleqpScalingData* scaling)
{
  SleqpProblem* problem = scaling->problem;
  SleqpProblem* scaled_problem = scaling->scaled_problem;

  SLEQP_CALL(sleqp_sparse_vector_copy(problem->var_lb,
                                      scaled_problem->var_lb));

  SLEQP_CALL(apply_unscaling(scaled_problem->var_lb,
                             scaling->var_weights,
                             0));

  SLEQP_CALL(sleqp_sparse_vector_copy(problem->var_ub,
                                      scaled_problem->var_ub));

  SLEQP_CALL(apply_unscaling(scaled_problem->var_ub,
                             scaling->var_weights,
                             0));

  SLEQP_CALL(sleqp_sparse_vector_copy(problem->cons_lb,
                                      scaled_problem->cons_lb));

  SLEQP_CALL(apply_scaling(scaled_problem->cons_lb,
                           scaling->cons_weights,
                           0));

  SLEQP_CALL(sleqp_sparse_vector_copy(problem->cons_ub,
                                      scaled_problem->cons_ub));

  SLEQP_CALL(apply_scaling(scaled_problem->cons_ub,
                           scaling->cons_weights,
                           0));

  return SLEQP_OKAY;
}

double sleqp_scale_func_val(SleqpScalingData* scaling,
                            double func_val)
{
  return ldexp(func_val, scaling->func_scale);
}

SLEQP_RETCODE sleqp_scale_point(SleqpScalingData* scaling,
                                SleqpSparseVec* point)
{
  SLEQP_CALL(apply_unscaling(point, scaling->var_weights, 0));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_scale_func_grad(SleqpScalingData* scaling,
                                    SleqpSparseVec* func_grad)
{
  SLEQP_CALL(apply_scaling(func_grad,
                           scaling->var_weights,
                           scaling->func_scale));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_scale_cons_val(SleqpScalingData* scaling,
                                   SleqpSparseVec* cons_val)
{
  SLEQP_CALL(apply_scaling(cons_val, scaling->cons_weights, 0));

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
                                  scaling->cons_weights[row] +
                                  scaling->var_weights[col]);
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_scale_cons_duals(SleqpScalingData* scaling,
                                     SleqpSparseVec* cons_duals)
{
  SLEQP_CALL(apply_unscaling(cons_duals,
                             scaling->cons_weights,
                             -1* scaling->func_scale));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_scale_var_duals(SleqpScalingData* scaling,
                                    SleqpSparseVec* var_duals)
{
  SLEQP_CALL(apply_scaling(var_duals,
                           scaling->var_weights,
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


double sleqp_unscale_func_val(SleqpScalingData* scaling,
                              double scaled_func_val)
{
  return ldexp(scaled_func_val, -1*scaling->func_scale);
}

SLEQP_RETCODE sleqp_unscale_point(SleqpScalingData* scaling,
                                  SleqpSparseVec* scaled_point)
{
  SLEQP_CALL(apply_scaling(scaled_point, scaling->var_weights, 0));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_unscale_func_grad(SleqpScalingData* scaling,
                                      SleqpSparseVec* scaled_func_grad)
{
  SLEQP_CALL(apply_unscaling(scaled_func_grad,
                             scaling->var_weights,
                             (-1)*scaling->func_scale));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_unscale_cons_val(SleqpScalingData* scaling,
                                     SleqpSparseVec* scaled_cons_val)
{
  SLEQP_CALL(apply_unscaling(scaled_cons_val, scaling->cons_weights, 0));

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
                                         - (scaling->var_weights[col] +
                                            scaling->cons_weights[row]));
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_unscale_cons_duals(SleqpScalingData* scaling,
                                       SleqpSparseVec* scaled_cons_duals)
{
  SLEQP_CALL(apply_scaling(scaled_cons_duals,
                           scaling->cons_weights,
                           -1* scaling->func_scale));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_unscale_var_duals(SleqpScalingData* scaling,
                                      SleqpSparseVec* scaled_var_duals)
{
  SLEQP_CALL(apply_unscaling(scaled_var_duals,
                             scaling->var_weights,
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

SLEQP_RETCODE sleqp_func_scaling_from_gradient(SleqpScalingData* scaling,
                                               SleqpSparseVec* gradient)
{
  double max_val = 0.;

  int* scaling_factor = &(scaling->func_scale);

  const double eps = sleqp_params_get_eps(scaling->params);

  *scaling_factor = 0;

  for(int k = 0; k < gradient->nnz; ++k)
  {
    double cur_val = SLEQP_ABS(gradient->data[k]);

    if(sleqp_zero(cur_val, eps))
    {
      continue;
    }

    max_val = SLEQP_MAX(max_val, cur_val);
  }

  if(max_val != 0.)
  {
    assert(max_val > 0.);

    frexp(1. / max_val, scaling_factor);

    --(*scaling_factor);

    CLIP_WEIGHT(*scaling_factor);

    assert((*scaling_factor) >= min_weight);
    assert((*scaling_factor) <= max_weight);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE compute_scale_factors(SleqpScalingData* scaling,
                                           SleqpSparseMatrix* matrix,
                                           bool column,
                                           int* scaling_factors,
                                           int* prescaling)
{
  double* max_vals = scaling->max_cache;

  const double eps = sleqp_params_get_eps(scaling->params);

  int size = (column) ? matrix->num_cols : matrix->num_rows;

  for(int i = 0; i < size; ++i)
  {
    max_vals[i] = 0.;

    scaling_factors[i] = 0;
  }

  int col = 0;

  for(int index = 0; index < matrix->nnz; ++index)
  {
    while(index >= matrix->cols[col + 1])
    {
      ++col;
    }

    int row = matrix->rows[index];

    int cur_entry = (column) ? col : row;
    int other_entry = (column) ? row : col;

    double cur_val = ldexp(matrix->data[index], prescaling[other_entry]);

    cur_val = SLEQP_ABS(cur_val);

    if(sleqp_zero(cur_val, eps))
    {
      continue;
    }

    max_vals[cur_entry] = SLEQP_MAX(max_vals[cur_entry], cur_val);
  }

  for(int i = 0; i < size; ++i)
  {
    if(max_vals[i] == 0.)
    {
      scaling_factors[i] = 0;
      continue;
    }

    assert(max_vals[i] > 0.);

    frexp(1.0 / (max_vals[i]), &(scaling_factors[i]));

    --scaling_factors[i];

    CLIP_WEIGHT(scaling_factors[i]);

    assert(scaling_factors[i] >= min_weight);
    assert(scaling_factors[i] <= max_weight);
  }

  return SLEQP_OKAY;
}

static double max_matrix_ratio(SleqpScalingData* scaling,
                               SleqpSparseMatrix* matrix,
                               bool column)
{
  double* min_vals = scaling->min_cache;
  double* max_vals = scaling->max_cache;

  const double eps = sleqp_params_get_eps(scaling->params);

  int size = (column) ? matrix->num_cols : matrix->num_rows;

  for(int i = 0; i < size; ++i)
  {
    min_vals[i] = sleqp_infinity();
    max_vals[i] = 0.;
  }

  int col = 0;

  for(int index = 0; index < matrix->nnz; ++index)
  {
    while(index >= matrix->cols[col + 1])
    {
      ++col;
    }

    int row = matrix->rows[index];
    double cur_val = SLEQP_ABS(matrix->data[index]);

    if(sleqp_zero(cur_val, eps))
    {
      continue;
    }

    int entry = (column) ? col : row;

    min_vals[entry] = SLEQP_MIN(min_vals[entry], cur_val);

    max_vals[entry] = SLEQP_MAX(max_vals[entry], cur_val);
  }

  double max_ratio = 0.;

  for(int i = 0; i < size; ++i)
  {
    if(max_vals[i] == 0.)
    {
      continue;
    }

    assert(min_vals[i] > 0.);
    assert(max_vals[i] > 0.);

    max_ratio = SLEQP_MAX(max_ratio, max_vals[i] / min_vals[i]);
  }

  return max_ratio;
}

SLEQP_RETCODE sleqp_scaling_from_cons_jac(SleqpScalingData* scaling,
                                          SleqpSparseMatrix* cons_jac)
{
  SLEQP_CALL(reset_scaling(scaling));

  SleqpProblem* problem = scaling->problem;

  if(!(scaling->min_cache))
  {
    const int size = SLEQP_MAX(problem->num_constraints,
                               problem->num_variables);

    SLEQP_CALL(sleqp_calloc(&scaling->min_cache, size));
    SLEQP_CALL(sleqp_calloc(&scaling->max_cache, size));
  }

  const double var_ratio = max_matrix_ratio(scaling, cons_jac, true);
  const double cons_ratio = max_matrix_ratio(scaling, cons_jac, false);

  const bool vars_first = var_ratio < cons_ratio;

  if(vars_first)
  {
    SLEQP_CALL(compute_scale_factors(scaling,
                                     cons_jac,
                                     true,
                                     scaling->var_weights,
                                     scaling->cons_weights));

    SLEQP_CALL(compute_scale_factors(scaling,
                                     cons_jac,
                                     false,
                                     scaling->cons_weights,
                                     scaling->var_weights));
  }
  else
  {
    SLEQP_CALL(compute_scale_factors(scaling,
                                     cons_jac,
                                     false,
                                     scaling->cons_weights,
                                     scaling->var_weights));

    SLEQP_CALL(compute_scale_factors(scaling,
                                     cons_jac,
                                     true,
                                     scaling->var_weights,
                                     scaling->cons_weights));
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_scaling_free(SleqpScalingData** star)
{
  SleqpScalingData* scaling = *star;

  if(!scaling)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_problem_free(&(scaling->scaled_problem)));

  SLEQP_CALL(sleqp_func_free(&(scaling->scaled_func)));

  SLEQP_CALL(sleqp_sparse_vector_free(&(scaling->scaled_cons_duals)));

  SLEQP_CALL(sleqp_sparse_vector_free(&(scaling->scaled_direction)));

  SLEQP_CALL(sleqp_sparse_vector_free(&(scaling->unscaled_value)));

  sleqp_free(&(scaling->max_cache));

  sleqp_free(&(scaling->min_cache));

  sleqp_free(&(scaling->cons_weights));

  sleqp_free(&(scaling->var_weights));

  sleqp_free(star);

  *star = NULL;

  return SLEQP_OKAY;
}
