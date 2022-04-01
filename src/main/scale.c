#include "scale.h"

#include <assert.h>
#include <fenv.h>
#include <math.h>

#include "cmp.h"
#include "error.h"
#include "log.h"
#include "mem.h"

#define MAX_WEIGHT 30

const int max_weight = MAX_WEIGHT;
const int min_weight = -(MAX_WEIGHT);

#define CLIP_WEIGHT(w) (w) = SLEQP_MIN(SLEQP_MAX((w), min_weight), max_weight)

struct SleqpScaling
{
  int refcount;

  int num_variables;
  int num_constraints;

  int obj_weight;
  int* var_weights;
  int* cons_weights;

  double* min_cache;
  double* max_cache;
};

static SLEQP_RETCODE
apply_const_scaling(SleqpVec* vec, int value)
{
  for (int k = 0; k < vec->nnz; ++k)
  {
    vec->data[k] = ldexp(vec->data[k], value);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
apply_const_unscaling(SleqpVec* vec, int value)
{
  SLEQP_CALL(apply_const_scaling(vec, (-1) * value));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
apply_unscaling(SleqpVec* vec, int* scales, int offset)
{
  for (int k = 0; k < vec->nnz; ++k)
  {
    vec->data[k] = ldexp(vec->data[k], scales[vec->indices[k]] + offset);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
apply_scaling(SleqpVec* vec, int* scales, int offset)
{
  for (int k = 0; k < vec->nnz; ++k)
  {
    vec->data[k] = ldexp(vec->data[k], -1 * scales[vec->indices[k]] + offset);
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_scaling_reset(SleqpScaling* scaling)
{
  for (int j = 0; j < scaling->num_variables; ++j)
  {
    scaling->var_weights[j] = 0;
  }

  for (int i = 0; i < scaling->num_constraints; ++i)
  {
    scaling->cons_weights[i] = 0;
  }

  scaling->obj_weight = 0;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_scaling_create(SleqpScaling** star,
                     int num_variables,
                     int num_constraints)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpScaling* scaling = *star;

  *scaling = (SleqpScaling){0};

  if (num_variables < 0)
  {
    sleqp_raise(SLEQP_ILLEGAL_ARGUMENT,
                "Negative number of variables (%d) provided to scaling",
                num_variables);
  }

  if (num_constraints < 0)
  {
    sleqp_raise(SLEQP_ILLEGAL_ARGUMENT,
                "Negative number of constraints (%d) provided to scaling",
                num_constraints);
  }

  scaling->refcount = 1;

  scaling->num_variables   = num_variables;
  scaling->num_constraints = num_constraints;

  SLEQP_CALL(
    sleqp_alloc_array(&(scaling->var_weights), scaling->num_variables));

  SLEQP_CALL(
    sleqp_alloc_array(&(scaling->cons_weights), scaling->num_constraints));

  scaling->obj_weight = 0;

  SLEQP_CALL(sleqp_scaling_reset(scaling));

  return SLEQP_OKAY;
}

int
sleqp_scaling_num_vars(SleqpScaling* scaling)
{
  return scaling->num_variables;
}

int
sleqp_scaling_num_cons(SleqpScaling* scaling)
{
  return scaling->num_constraints;
}

int
sleqp_scaling_obj_weight(SleqpScaling* scaling)
{
  return scaling->obj_weight;
}

SLEQP_RETCODE
sleqp_scaling_set_obj_weight(SleqpScaling* scaling, int weight)
{
  scaling->obj_weight = weight;
  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_scaling_set_obj_weight_from_nominal(SleqpScaling* scaling,
                                          double nominal_value)
{
  frexp(nominal_value, &(scaling->obj_weight));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_scaling_set_var_weight(SleqpScaling* scaling, int index, int weight)
{
  if (index < 0 || index >= scaling->num_variables)
  {
    sleqp_raise(SLEQP_ILLEGAL_ARGUMENT, "Invalid variable index %d", index);
  }

  scaling->var_weights[index] = weight;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_scaling_set_var_weights_from_nominal(SleqpScaling* scaling,
                                           double* nominal_values)
{
  for (int j = 0; j < scaling->num_variables; ++j)
  {
    frexp(nominal_values[j], scaling->var_weights + j);
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_scaling_set_var_weight_from_nominal(SleqpScaling* scaling,
                                          int index,
                                          double nominal_value)
{
  if (index < 0 || index >= scaling->num_variables)
  {
    sleqp_raise(SLEQP_ILLEGAL_ARGUMENT, "Invalid variable index %d", index);
  }

  frexp(nominal_value, scaling->var_weights + index);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_scaling_set_cons_weight(SleqpScaling* scaling, int index, int weight)
{
  if (index < 0 || index >= scaling->num_constraints)
  {
    sleqp_raise(SLEQP_ILLEGAL_ARGUMENT, "Invalid constraint index %d", index);
  }

  scaling->cons_weights[index] = weight;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_scaling_set_cons_weights_from_nominal(SleqpScaling* scaling,
                                            double* nominal_values)
{
  for (int i = 0; i < scaling->num_constraints; ++i)
  {
    frexp(nominal_values[i], scaling->cons_weights + i);
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_scaling_set_cons_weight_from_nominal(SleqpScaling* scaling,
                                           int index,
                                           double nominal_value)
{
  if (index < 0 || index >= scaling->num_constraints)
  {
    sleqp_raise(SLEQP_ILLEGAL_ARGUMENT, "Invalid constraint index %d", index);
  }

  frexp(nominal_value, scaling->cons_weights + index);

  return SLEQP_OKAY;
}

int*
sleqp_scaling_var_weights(SleqpScaling* scaling)
{
  return scaling->var_weights;
}

int*
sleqp_scaling_cons_weights(SleqpScaling* scaling)
{
  return scaling->cons_weights;
}

double
sleqp_scale_obj_val(SleqpScaling* scaling, double obj_val)
{
  return ldexp(obj_val, (-1) * scaling->obj_weight);
}

double
sleqp_scale_lsq_obj_val(SleqpScaling* scaling, double obj_val)
{
  return ldexp(obj_val, (-1) * 2 * scaling->obj_weight);
}

SLEQP_RETCODE
sleqp_scale_lsq_residuals(SleqpScaling* scaling, SleqpVec* lsq_residuals)
{
  return apply_const_unscaling(lsq_residuals, scaling->obj_weight);
}

SLEQP_RETCODE
sleqp_scale_lsq_forward_direction(SleqpScaling* scaling,
                                  SleqpVec* forward_direction)
{
  return apply_unscaling(forward_direction, scaling->var_weights, 0);
}

SLEQP_RETCODE
sleqp_scale_lsq_adjoint_direction(SleqpScaling* scaling,
                                  SleqpVec* adjoint_direction)
{
  return apply_const_unscaling(adjoint_direction, scaling->obj_weight);
}

SLEQP_RETCODE
sleqp_scale_point(SleqpScaling* scaling, SleqpVec* point)
{
  SLEQP_CALL(apply_scaling(point, scaling->var_weights, 0));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_scale_obj_grad(SleqpScaling* scaling, SleqpVec* obj_grad)
{
  SLEQP_CALL(apply_unscaling(obj_grad,
                             scaling->var_weights,
                             (-1) * scaling->obj_weight));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_scale_cons_val(SleqpScaling* scaling, SleqpVec* cons_val)
{
  SLEQP_CALL(apply_scaling(cons_val, scaling->cons_weights, 0));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_scale_cons_general(SleqpScaling* scaling, SleqpVec* general_val)
{
  SLEQP_CALL(apply_scaling(general_val, scaling->cons_weights, 0));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_scale_cons_linear(SleqpScaling* scaling, SleqpVec* linear_val)
{
  const int num_linear  = linear_val->dim;
  const int num_cons    = scaling->num_constraints;
  const int num_general = num_cons - num_linear;

  SLEQP_CALL(apply_scaling(linear_val, scaling->cons_weights + num_general, 0));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_scale_cons_jac(SleqpScaling* scaling, SleqpSparseMatrix* cons_jac)
{
  int col = 0;

  const int* cons_jac_cols = sleqp_sparse_matrix_cols(cons_jac);
  const int* cons_jac_rows = sleqp_sparse_matrix_rows(cons_jac);
  double* cons_jac_data    = sleqp_sparse_matrix_data(cons_jac);

  int cons_jac_nnz = sleqp_sparse_matrix_nnz(cons_jac);

  for (int index = 0; index < cons_jac_nnz; ++index)
  {
    while (index >= cons_jac_cols[col + 1])
    {
      ++col;
    }

    const int row = cons_jac_rows[index];

    cons_jac_data[index]
      = ldexp(cons_jac_data[index],
              (-1) * scaling->cons_weights[row] + scaling->var_weights[col]);
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_scale_linear_coeffs(SleqpScaling* scaling,
                          SleqpSparseMatrix* linear_coeffs)
{
  int col = 0;

  const int num_linear = sleqp_sparse_matrix_num_rows(linear_coeffs);

  assert(num_linear <= scaling->num_constraints);

  const int* linear_coeffs_cols = sleqp_sparse_matrix_cols(linear_coeffs);
  const int* linear_coeffs_rows = sleqp_sparse_matrix_rows(linear_coeffs);
  double* linear_coeffs_data    = sleqp_sparse_matrix_data(linear_coeffs);

  int linear_coeffs_nnz = sleqp_sparse_matrix_nnz(linear_coeffs);

  for (int index = 0; index < linear_coeffs_nnz; ++index)
  {
    while (index >= linear_coeffs_cols[col + 1])
    {
      ++col;
    }

    const int row = linear_coeffs_rows[index];

    linear_coeffs_data[index]
      = ldexp(linear_coeffs_data[index],
              (-1) * scaling->cons_weights[num_linear + row]
                + scaling->var_weights[col]);
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_scale_cons_duals(SleqpScaling* scaling, SleqpVec* cons_duals)
{
  SLEQP_CALL(apply_unscaling(cons_duals,
                             scaling->cons_weights,
                             (-1) * scaling->obj_weight));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_scale_var_duals(SleqpScaling* scaling, SleqpVec* var_duals)
{
  SLEQP_CALL(apply_unscaling(var_duals,
                             scaling->var_weights,
                             (-1) * scaling->obj_weight));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_scale_hessian_product(SleqpScaling* scaling, SleqpVec* product)
{
  SLEQP_CALL(
    apply_unscaling(product, scaling->var_weights, (-1) * scaling->obj_weight));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_scale_iterate(SleqpScaling* scaling,
                    SleqpIterate* scaled_iterate,
                    bool lsq)
{
  SLEQP_CALL(sleqp_scale_point(scaling, sleqp_iterate_primal(scaled_iterate)));

  const double obj_val = sleqp_iterate_obj_val(scaled_iterate);

  if (lsq)
  {
    const double scaled_obj_val = sleqp_scale_lsq_obj_val(scaling, obj_val);

    SLEQP_CALL(sleqp_iterate_set_obj_val(scaled_iterate, scaled_obj_val));
  }
  else
  {
    const double scaled_obj_val = sleqp_scale_obj_val(scaling, obj_val);

    SLEQP_CALL(sleqp_iterate_set_obj_val(scaled_iterate, scaled_obj_val));
  }

  SLEQP_CALL(
    sleqp_scale_obj_grad(scaling, sleqp_iterate_obj_grad(scaled_iterate)));

  SLEQP_CALL(
    sleqp_scale_cons_val(scaling, sleqp_iterate_cons_val(scaled_iterate)));

  SLEQP_CALL(
    sleqp_scale_cons_jac(scaling, sleqp_iterate_cons_jac(scaled_iterate)));

  SLEQP_CALL(
    sleqp_scale_cons_duals(scaling, sleqp_iterate_cons_dual(scaled_iterate)));

  SLEQP_CALL(
    sleqp_scale_var_duals(scaling, sleqp_iterate_vars_dual(scaled_iterate)));

  return SLEQP_OKAY;
}

double
sleqp_unscale_obj_val(SleqpScaling* scaling, double scaled_obj_val)
{
  return ldexp(scaled_obj_val, scaling->obj_weight);
}

double
sleqp_unscale_lsq_obj_val(SleqpScaling* scaling, double scaled_obj_val)
{
  return ldexp(scaled_obj_val, 2 * scaling->obj_weight);
}

SLEQP_RETCODE
sleqp_unscale_point(SleqpScaling* scaling, SleqpVec* scaled_point)
{
  SLEQP_CALL(apply_unscaling(scaled_point, scaling->var_weights, 0));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_unscale_obj_grad(SleqpScaling* scaling, SleqpVec* scaled_obj_grad)
{
  SLEQP_CALL(
    apply_scaling(scaled_obj_grad, scaling->var_weights, scaling->obj_weight));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_unscale_cons_val(SleqpScaling* scaling, SleqpVec* scaled_cons_val)
{
  SLEQP_CALL(apply_unscaling(scaled_cons_val, scaling->cons_weights, 0));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_unscale_cons_jac(SleqpScaling* scaling,
                       SleqpSparseMatrix* scaled_cons_jac)
{
  int col = 0;

  const int* cons_jac_cols = sleqp_sparse_matrix_cols(scaled_cons_jac);
  const int* cons_jac_rows = sleqp_sparse_matrix_rows(scaled_cons_jac);
  double* cons_jac_data    = sleqp_sparse_matrix_data(scaled_cons_jac);

  int cons_jac_nnz = sleqp_sparse_matrix_nnz(scaled_cons_jac);

  for (int index = 0; index < cons_jac_nnz; ++index)
  {
    while (index >= cons_jac_cols[col + 1])
    {
      ++col;
    }

    const int row = cons_jac_rows[index];

    cons_jac_data[index]
      = ldexp(cons_jac_data[index],
              (-1) * (scaling->var_weights[col]) + scaling->cons_weights[row]);
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_unscale_cons_duals(SleqpScaling* scaling, SleqpVec* scaled_cons_duals)
{
  SLEQP_CALL(apply_scaling(scaled_cons_duals,
                           scaling->cons_weights,
                           scaling->obj_weight));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_unscale_var_duals(SleqpScaling* scaling, SleqpVec* scaled_var_duals)
{
  SLEQP_CALL(
    apply_scaling(scaled_var_duals, scaling->var_weights, scaling->obj_weight));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_unscale_hessian_direction(SleqpScaling* scaling,
                                SleqpVec* direction,
                                SleqpVec* cons_duals)
{
  SLEQP_CALL(apply_unscaling(direction, scaling->var_weights, 0));

  SLEQP_CALL(
    apply_scaling(cons_duals, scaling->cons_weights, scaling->obj_weight));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_unscale_iterate(SleqpScaling* scaling,
                      SleqpIterate* scaled_iterate,
                      bool lsq)
{
  SLEQP_CALL(
    sleqp_unscale_point(scaling, sleqp_iterate_primal(scaled_iterate)));

  const double obj_val = sleqp_iterate_obj_val(scaled_iterate);

  if (lsq)
  {
    const double unscaled_obj_val = sleqp_unscale_lsq_obj_val(scaling, obj_val);

    SLEQP_CALL(sleqp_iterate_set_obj_val(scaled_iterate, unscaled_obj_val));
  }
  else
  {
    const double unscaled_obj_val = sleqp_unscale_obj_val(scaling, obj_val);

    SLEQP_CALL(sleqp_iterate_set_obj_val(scaled_iterate, unscaled_obj_val));
  }

  SLEQP_CALL(
    sleqp_unscale_obj_grad(scaling, sleqp_iterate_obj_grad(scaled_iterate)));

  SLEQP_CALL(
    sleqp_unscale_cons_val(scaling, sleqp_iterate_cons_val(scaled_iterate)));

  SLEQP_CALL(
    sleqp_unscale_cons_jac(scaling, sleqp_iterate_cons_jac(scaled_iterate)));

  SLEQP_CALL(
    sleqp_unscale_cons_duals(scaling, sleqp_iterate_cons_dual(scaled_iterate)));

  SLEQP_CALL(
    sleqp_unscale_var_duals(scaling, sleqp_iterate_vars_dual(scaled_iterate)));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_obj_scaling_from_grad(SleqpScaling* scaling,
                            SleqpVec* gradient,
                            double eps)
{
  double max_val = 0.;

  int* scaling_factor = &(scaling->obj_weight);

  *scaling_factor = 0;

  for (int k = 0; k < gradient->nnz; ++k)
  {
    double cur_val = SLEQP_ABS(gradient->data[k]);

    if (sleqp_is_zero(cur_val, eps))
    {
      continue;
    }

    max_val = SLEQP_MAX(max_val, cur_val);
  }

  if (max_val != 0.)
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

static SLEQP_RETCODE
compute_scale_factors(SleqpScaling* scaling,
                      SleqpSparseMatrix* matrix,
                      bool column,
                      int* scaling_factors,
                      int* prescaling,
                      double eps)
{
  double* max_vals = scaling->max_cache;

  const int num_cols = sleqp_sparse_matrix_num_cols(matrix);
  const int num_rows = sleqp_sparse_matrix_num_rows(matrix);

  int size = (column) ? num_cols : num_rows;

  for (int i = 0; i < size; ++i)
  {
    max_vals[i] = 0.;

    scaling_factors[i] = 0;
  }

  int col = 0;

  const int* matrix_cols = sleqp_sparse_matrix_cols(matrix);
  const int* matrix_rows = sleqp_sparse_matrix_rows(matrix);
  double* matrix_data    = sleqp_sparse_matrix_data(matrix);

  int matrix_nnz = sleqp_sparse_matrix_nnz(matrix);

  for (int index = 0; index < matrix_nnz; ++index)
  {
    while (index >= matrix_cols[col + 1])
    {
      ++col;
    }

    int row = matrix_rows[index];

    int cur_entry   = (column) ? col : row;
    int other_entry = (column) ? row : col;

    double cur_val = ldexp(matrix_data[index], prescaling[other_entry]);

    cur_val = SLEQP_ABS(cur_val);

    if (sleqp_is_zero(cur_val, eps))
    {
      continue;
    }

    max_vals[cur_entry] = SLEQP_MAX(max_vals[cur_entry], cur_val);
  }

  for (int i = 0; i < size; ++i)
  {
    if (max_vals[i] == 0.)
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

static double
max_matrix_ratio(SleqpScaling* scaling,
                 SleqpSparseMatrix* matrix,
                 bool column,
                 double eps)
{
  double* min_vals = scaling->min_cache;
  double* max_vals = scaling->max_cache;

  const int num_cols = sleqp_sparse_matrix_num_cols(matrix);
  const int num_rows = sleqp_sparse_matrix_num_rows(matrix);

  const int* matrix_cols = sleqp_sparse_matrix_cols(matrix);
  const int* matrix_rows = sleqp_sparse_matrix_rows(matrix);
  double* matrix_data    = sleqp_sparse_matrix_data(matrix);

  const int matrix_nnz = sleqp_sparse_matrix_nnz(matrix);

  int size = (column) ? num_cols : num_rows;

  for (int i = 0; i < size; ++i)
  {
    min_vals[i] = sleqp_infinity();
    max_vals[i] = 0.;
  }

  int col = 0;

  for (int index = 0; index < matrix_nnz; ++index)
  {
    while (index >= matrix_cols[col + 1])
    {
      ++col;
    }

    int row        = matrix_rows[index];
    double cur_val = SLEQP_ABS(matrix_data[index]);

    if (sleqp_is_zero(cur_val, eps))
    {
      continue;
    }

    int entry = (column) ? col : row;

    min_vals[entry] = SLEQP_MIN(min_vals[entry], cur_val);

    max_vals[entry] = SLEQP_MAX(max_vals[entry], cur_val);
  }

  double max_ratio = 0.;

  for (int i = 0; i < size; ++i)
  {
    if (max_vals[i] == 0.)
    {
      continue;
    }

    assert(min_vals[i] > 0.);
    assert(max_vals[i] > 0.);

    max_ratio = SLEQP_MAX(max_ratio, max_vals[i] / min_vals[i]);
  }

  return max_ratio;
}

SLEQP_RETCODE
sleqp_scaling_from_cons_jac(SleqpScaling* scaling,
                            SleqpSparseMatrix* cons_jac,
                            double eps)
{
  SLEQP_CALL(sleqp_scaling_reset(scaling));

  if (!(scaling->min_cache))
  {
    const int size
      = SLEQP_MAX(scaling->num_constraints, scaling->num_variables);

    SLEQP_CALL(sleqp_alloc_array(&scaling->min_cache, size));
    SLEQP_CALL(sleqp_alloc_array(&scaling->max_cache, size));
  }

  const double var_ratio  = max_matrix_ratio(scaling, cons_jac, true, eps);
  const double cons_ratio = max_matrix_ratio(scaling, cons_jac, false, eps);

  const bool vars_first = var_ratio < cons_ratio;

  if (vars_first)
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

static SLEQP_RETCODE
scaling_free(SleqpScaling** star)
{
  SleqpScaling* scaling = *star;

  if (!scaling)
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

SLEQP_RETCODE
sleqp_scaling_capture(SleqpScaling* scaling)
{
  ++scaling->refcount;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_scaling_release(SleqpScaling** star)
{
  SleqpScaling* scaling = *star;

  if (!scaling)
  {
    return SLEQP_OKAY;
  }

  if (--scaling->refcount == 0)
  {
    SLEQP_CALL(scaling_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
