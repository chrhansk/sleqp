#include "sleqp_scale.h"

#include <math.h>
#include <fenv.h>

#include "sleqp_cmp.h"
#include "sleqp_mem.h"

#define INIT_MATH_ERROR_CHECK                   \
  fenv_t fenv_current;                          \
  if(math_errhandling & MATH_ERREXCEPT)         \
  {                                             \
    fegetenv(&fenv_current);                    \
  }                                             \
  while(false)

#define MATH_ERROR_CHECK(error_flags)                  \
  if(math_errhandling & MATH_ERREXCEPT)                \
  {                                                    \
    const bool has_error = fetestexcept(error_flags);  \
                                                       \
    fesetenv(&fenv_current);                           \
                                                       \
    if(has_error)                                      \
    {                                                  \
      return SLEQP_MATH_ERROR;                         \
    }                                                  \
  }                                                    \
  while(false)

#define SCALING_ERROR_FLAGS (FE_OVERFLOW | FE_UNDERFLOW)


const int max_weight = 30;
const int min_weight = -max_weight;

#define CLIP_WEIGHT(w) \
  (w) = SLEQP_MIN(SLEQP_MAX((w), min_weight), max_weight)

struct SleqpScalingData
{
  int refcount;

  int num_variables;
  int num_constraints;

  int func_weight;
  int* var_weights;
  int* cons_weights;

  double* min_cache;
  double* max_cache;
};

static SLEQP_RETCODE apply_scaling(SleqpSparseVec* vec,
                                   int* scales,
                                   int offset)
{
  INIT_MATH_ERROR_CHECK;

  for(int k = 0; k < vec->nnz; ++k)
  {
    vec->data[k] = ldexp(vec->data[k], scales[vec->indices[k]] + offset);
  }

  MATH_ERROR_CHECK(SCALING_ERROR_FLAGS);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE apply_unscaling(SleqpSparseVec* vec,
                                     int* scales,
                                     int offset)
{
  INIT_MATH_ERROR_CHECK;

  for(int k = 0; k < vec->nnz; ++k)
  {
    vec->data[k] = ldexp(vec->data[k], -1*scales[vec->indices[k]] + offset);
  }

  MATH_ERROR_CHECK(SCALING_ERROR_FLAGS);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE reset_scaling(SleqpScalingData* scaling)
{
  for(int j = 0; j < scaling->num_variables; ++j)
  {
    scaling->var_weights[j] = 0;
  }

  for(int i = 0; i < scaling->num_constraints; ++i)
  {
    scaling->cons_weights[i] = 0;
  }

  scaling->func_weight = 0;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_scaling_create(SleqpScalingData** star,
                                   int num_variables,
                                   int num_constraints)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpScalingData* scaling = *star;

  *scaling = (SleqpScalingData) {0};

  if(num_variables < 0)
  {
    sleqp_log_error("Negative number of variables provided to scaling");
    return SLEQP_ILLEGAL_ARGUMENT;
  }

  if(num_constraints < 0)
  {
    sleqp_log_error("Negative number of constraints provided to scaling");
    return SLEQP_ILLEGAL_ARGUMENT;
  }

  scaling->refcount = 1;

  scaling->num_variables = num_variables;
  scaling->num_constraints = num_constraints;

  SLEQP_CALL(sleqp_calloc(&(scaling->var_weights),
                          scaling->num_variables));

  SLEQP_CALL(sleqp_calloc(&(scaling->cons_weights),
                          scaling->num_constraints));

  scaling->func_weight = 0;

  SLEQP_CALL(reset_scaling(scaling));

  return SLEQP_OKAY;
}

int sleqp_scaling_get_num_variables(SleqpScalingData* scaling)
{
  return scaling->num_variables;
}

int sleqp_scaling_get_num_constraints(SleqpScalingData* scaling)
{
  return scaling->num_constraints;
}

int sleqp_scaling_get_func_weight(SleqpScalingData* scaling)
{
  return scaling->func_weight;
}

SLEQP_RETCODE sleqp_scaling_set_func_weight(SleqpScalingData* scaling,
                                            int weight)
{
  scaling->func_weight = weight;
  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_scaling_set_func_weight_from_nominal(SleqpScalingData* scaling,
                                                         double nominal_value)
{
  frexp(nominal_value, &(scaling->func_weight));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_scaling_set_var_weight(SleqpScalingData* scaling,
                                           int index,
                                           int weight)
{
  if(index < 0 || index >= scaling->num_variables)
  {
    return SLEQP_ILLEGAL_ARGUMENT;
  }

  scaling->var_weights[index] = weight;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_scaling_set_var_weights_from_nominal(SleqpScalingData* scaling,
                                                         double* nominal_values)
{
  for(int j = 0; j < scaling->num_variables; ++j)
  {
    frexp(nominal_values[j], scaling->var_weights + j);
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_scaling_set_var_weight_from_nominal(SleqpScalingData* scaling,
                                                        int index,
                                                        double nominal_value)
{
  if(index < 0 || index >= scaling->num_variables)
  {
    return SLEQP_ILLEGAL_ARGUMENT;
  }

  frexp(nominal_value, scaling->var_weights + index);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_scaling_set_cons_weight(SleqpScalingData* scaling,
                                            int index,
                                            int weight)
{
  if(index < 0 || index >= scaling->num_constraints)
  {
    return SLEQP_ILLEGAL_ARGUMENT;
  }

  scaling->cons_weights[index] = weight;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_scaling_set_cons_weights_from_nominal(SleqpScalingData* scaling,
                                                          double* nominal_values)
{
  for(int i = 0; i < scaling->num_constraints; ++i)
  {
    frexp(nominal_values[i], scaling->cons_weights + i);
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_scaling_set_cons_weight_from_nominal(SleqpScalingData* scaling,
                                                         int index,
                                                         double nominal_value)
{
  if(index < 0 || index >= scaling->num_constraints)
  {
    return SLEQP_ILLEGAL_ARGUMENT;
  }

  frexp(nominal_value, scaling->cons_weights + index);

  return SLEQP_OKAY;
}

int* sleqp_scaling_get_var_weights(SleqpScalingData* scaling)
{
  return scaling->var_weights;
}

int* sleqp_scaling_get_cons_weights(SleqpScalingData* scaling)
{
  return scaling->cons_weights;
}

double sleqp_scale_func_val(SleqpScalingData* scaling,
                            double func_val)
{
  return ldexp(func_val, (-1) * scaling->func_weight);
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
                           (-1) * scaling->func_weight));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_scale_cons_val(SleqpScalingData* scaling,
                                   SleqpSparseVec* cons_val)
{
  SLEQP_CALL(apply_unscaling(cons_val, scaling->cons_weights, 0));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_scale_cons_jac(SleqpScalingData* scaling,
                                   SleqpSparseMatrix* cons_jac)
{
  int col = 0;

  const int* cons_jac_cols = sleqp_sparse_matrix_get_cols(cons_jac);
  const int* cons_jac_rows = sleqp_sparse_matrix_get_rows(cons_jac);
  double* cons_jac_data = sleqp_sparse_matrix_get_data(cons_jac);

  int cons_jac_nnz = sleqp_sparse_matrix_get_nnz(cons_jac);

  INIT_MATH_ERROR_CHECK;

  for(int index = 0; index < cons_jac_nnz; ++index)
  {
    while(index >= cons_jac_cols[col + 1])
    {
      ++col;
    }

    const int row = cons_jac_rows[index];

    cons_jac_data[index] = ldexp(cons_jac_data[index],
                                 (-1) * scaling->cons_weights[row] +
                                 scaling->var_weights[col]);
  }

  MATH_ERROR_CHECK(SCALING_ERROR_FLAGS);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_scale_cons_duals(SleqpScalingData* scaling,
                                     SleqpSparseVec* cons_duals)
{
  SLEQP_CALL(apply_scaling(cons_duals,
                           scaling->cons_weights,
                           scaling->func_weight));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_scale_var_duals(SleqpScalingData* scaling,
                                    SleqpSparseVec* var_duals)
{
  SLEQP_CALL(apply_scaling(var_duals,
                           scaling->var_weights,
                           scaling->func_weight));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_scale_hessian_product(SleqpScalingData* scaling,
                                          SleqpSparseVec* product)
{
  SLEQP_CALL(apply_scaling(product,
                           scaling->var_weights,
                           (-1) * scaling->func_weight));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_scale_iterate(SleqpScalingData* scaling,
                                  SleqpIterate* scaled_iterate)
{
  SLEQP_CALL(sleqp_scale_point(scaling, sleqp_iterate_get_primal(scaled_iterate)));

  sleqp_iterate_set_func_val(scaled_iterate,
                             sleqp_scale_func_val(scaling, sleqp_iterate_get_func_val(scaled_iterate)));

  SLEQP_CALL(sleqp_scale_func_grad(scaling,
                                   sleqp_iterate_get_func_grad(scaled_iterate)));

  SLEQP_CALL(sleqp_scale_cons_val(scaling,
                                  sleqp_iterate_get_cons_val(scaled_iterate)));

  SLEQP_CALL(sleqp_scale_cons_jac(scaling,
                                  sleqp_iterate_get_cons_jac(scaled_iterate)));

  SLEQP_CALL(sleqp_scale_cons_duals(scaling,
                                    sleqp_iterate_get_cons_dual(scaled_iterate)));

  SLEQP_CALL(sleqp_scale_var_duals(scaling,
                                   sleqp_iterate_get_vars_dual(scaled_iterate)));

  return SLEQP_OKAY;
}


double sleqp_unscale_func_val(SleqpScalingData* scaling,
                              double scaled_func_val)
{
  return ldexp(scaled_func_val, scaling->func_weight);
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
                             scaling->func_weight));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_unscale_cons_val(SleqpScalingData* scaling,
                                     SleqpSparseVec* scaled_cons_val)
{
  SLEQP_CALL(apply_scaling(scaled_cons_val, scaling->cons_weights, 0));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_unscale_cons_jac(SleqpScalingData* scaling,
                                     SleqpSparseMatrix* scaled_cons_jac)
{
  int col = 0;

  const int* cons_jac_cols = sleqp_sparse_matrix_get_cols(scaled_cons_jac);
  const int* cons_jac_rows = sleqp_sparse_matrix_get_rows(scaled_cons_jac);
  double* cons_jac_data = sleqp_sparse_matrix_get_data(scaled_cons_jac);

  int cons_jac_nnz = sleqp_sparse_matrix_get_nnz(scaled_cons_jac);

  INIT_MATH_ERROR_CHECK;

  for(int index = 0; index < cons_jac_nnz; ++index)
  {
    while(index >= cons_jac_cols[col + 1])
    {
      ++col;
    }

    const int row = cons_jac_rows[index];

    cons_jac_data[index] = ldexp(cons_jac_data[index],
                                 (-1)*(scaling->var_weights[col]) +
                                 scaling->cons_weights[row]);
  }

  MATH_ERROR_CHECK(SCALING_ERROR_FLAGS);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_unscale_cons_duals(SleqpScalingData* scaling,
                                       SleqpSparseVec* scaled_cons_duals)
{
  SLEQP_CALL(apply_unscaling(scaled_cons_duals,
                             scaling->cons_weights,
                             (-1) * scaling->func_weight));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_unscale_var_duals(SleqpScalingData* scaling,
                                      SleqpSparseVec* scaled_var_duals)
{
  SLEQP_CALL(apply_unscaling(scaled_var_duals,
                             scaling->var_weights,
                             scaling->func_weight));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_unscale_hessian_direction(SleqpScalingData* scaling,
                                              SleqpSparseVec* direction,
                                              SleqpSparseVec* cons_duals)
{
  SLEQP_CALL(apply_scaling(direction,
                           scaling->var_weights,
                           0));

  SLEQP_CALL(apply_unscaling(cons_duals,
                             scaling->cons_weights,
                             scaling->func_weight));

  return SLEQP_OKAY;
}


SLEQP_RETCODE sleqp_unscale_iterate(SleqpScalingData* scaling,
                                    SleqpIterate* scaled_iterate)
{
  SLEQP_CALL(sleqp_unscale_point(scaling,
                                 sleqp_iterate_get_primal(scaled_iterate)));

  double func_val = sleqp_iterate_get_func_val(scaled_iterate);

  sleqp_iterate_set_func_val(scaled_iterate,
                             sleqp_unscale_func_val(scaling, func_val));

  SLEQP_CALL(sleqp_unscale_func_grad(scaling,
                                     sleqp_iterate_get_func_grad(scaled_iterate)));

  SLEQP_CALL(sleqp_unscale_cons_val(scaling,
                                    sleqp_iterate_get_cons_val(scaled_iterate)));

  SLEQP_CALL(sleqp_unscale_cons_jac(scaling,
                                    sleqp_iterate_get_cons_jac(scaled_iterate)));

  SLEQP_CALL(sleqp_unscale_cons_duals(scaling,
                                      sleqp_iterate_get_cons_dual(scaled_iterate)));

  SLEQP_CALL(sleqp_unscale_var_duals(scaling,
                                     sleqp_iterate_get_vars_dual(scaled_iterate)));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_func_scaling_from_gradient(SleqpScalingData* scaling,
                                               SleqpSparseVec* gradient,
                                               double eps)
{
  double max_val = 0.;

  int* scaling_factor = &(scaling->func_weight);

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
                                           int* prescaling,
                                           double eps)
{
  double* max_vals = scaling->max_cache;

  const int num_cols = sleqp_sparse_matrix_get_num_cols(matrix);
  const int num_rows = sleqp_sparse_matrix_get_num_rows(matrix);

  int size = (column) ? num_cols : num_rows;

  for(int i = 0; i < size; ++i)
  {
    max_vals[i] = 0.;

    scaling_factors[i] = 0;
  }

  int col = 0;

  const int* matrix_cols = sleqp_sparse_matrix_get_cols(matrix);
  const int* matrix_rows = sleqp_sparse_matrix_get_rows(matrix);
  double* matrix_data = sleqp_sparse_matrix_get_data(matrix);

  int matrix_nnz = sleqp_sparse_matrix_get_nnz(matrix);

  INIT_MATH_ERROR_CHECK;

  for(int index = 0; index < matrix_nnz; ++index)
  {
    while(index >= matrix_cols[col + 1])
    {
      ++col;
    }

    int row = matrix_rows[index];

    int cur_entry = (column) ? col : row;
    int other_entry = (column) ? row : col;

    double cur_val = ldexp(matrix_data[index], prescaling[other_entry]);

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

  MATH_ERROR_CHECK(SCALING_ERROR_FLAGS);

  return SLEQP_OKAY;
}

static double max_matrix_ratio(SleqpScalingData* scaling,
                               SleqpSparseMatrix* matrix,
                               bool column,
                               double eps)
{
  double* min_vals = scaling->min_cache;
  double* max_vals = scaling->max_cache;

  const int num_cols = sleqp_sparse_matrix_get_num_cols(matrix);
  const int num_rows = sleqp_sparse_matrix_get_num_rows(matrix);

  const int* matrix_cols = sleqp_sparse_matrix_get_cols(matrix);
  const int* matrix_rows = sleqp_sparse_matrix_get_rows(matrix);
  double* matrix_data = sleqp_sparse_matrix_get_data(matrix);

  const int matrix_nnz = sleqp_sparse_matrix_get_nnz(matrix);

  int size = (column) ? num_cols : num_rows;

  for(int i = 0; i < size; ++i)
  {
    min_vals[i] = sleqp_infinity();
    max_vals[i] = 0.;
  }

  int col = 0;

  for(int index = 0; index < matrix_nnz; ++index)
  {
    while(index >= matrix_cols[col + 1])
    {
      ++col;
    }

    int row = matrix_rows[index];
    double cur_val = SLEQP_ABS(matrix_data[index]);

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
                                          SleqpSparseMatrix* cons_jac,
                                          double eps)
{
  SLEQP_CALL(reset_scaling(scaling));

  if(!(scaling->min_cache))
  {
    const int size = SLEQP_MAX(scaling->num_constraints,
                               scaling->num_variables);

    SLEQP_CALL(sleqp_calloc(&scaling->min_cache, size));
    SLEQP_CALL(sleqp_calloc(&scaling->max_cache, size));
  }

  const double var_ratio = max_matrix_ratio(scaling, cons_jac, true, eps);
  const double cons_ratio = max_matrix_ratio(scaling, cons_jac, false, eps);

  const bool vars_first = var_ratio < cons_ratio;

  if(vars_first)
  {
    SLEQP_CALL(compute_scale_factors(scaling,
                                     cons_jac,
                                     true,
                                     scaling->var_weights,
                                     scaling->cons_weights,
                                     eps));

    SLEQP_CALL(compute_scale_factors(scaling,
                                     cons_jac,
                                     false,
                                     scaling->cons_weights,
                                     scaling->var_weights,
                                     eps));
  }
  else
  {
    SLEQP_CALL(compute_scale_factors(scaling,
                                     cons_jac,
                                     false,
                                     scaling->cons_weights,
                                     scaling->var_weights,
                                     eps));

    SLEQP_CALL(compute_scale_factors(scaling,
                                     cons_jac,
                                     true,
                                     scaling->var_weights,
                                     scaling->cons_weights,
                                     eps));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE scaling_free(SleqpScalingData** star)
{
  SleqpScalingData* scaling = *star;

  if(!scaling)
  {
    return SLEQP_OKAY;
  }

  sleqp_free(&(scaling->max_cache));

  sleqp_free(&(scaling->min_cache));

  sleqp_free(&(scaling->cons_weights));

  sleqp_free(&(scaling->var_weights));

  sleqp_free(star);

  *star = NULL;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_scaling_capture(SleqpScalingData* scaling)
{
  ++scaling->refcount;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_scaling_release(SleqpScalingData** star)
{
  SleqpScalingData* scaling = *star;

  if(!scaling)
  {
    return SLEQP_OKAY;
  }

  if(--scaling->refcount == 0)
  {
    SLEQP_CALL(scaling_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
