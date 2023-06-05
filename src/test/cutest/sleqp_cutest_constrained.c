#include "sleqp_cutest_constrained.h"

#include <stdlib.h>

#include "fail.h"
#include "log.h"
#include "mem.h"
#include "sparse/mat.h"

#include "sleqp_cutest_types.h"

typedef struct CUTestConsFuncData
{
  double eps;
  double zero_eps;

  int num_variables;
  int num_constraints;
  int num_linear;

  double* x;
  double* cons_vals;
  double* obj_grad;

  double* direction;
  double* multipliers;
  double* hessian_product;
  double* cons_hessian_product;

  int jac_nnz;
  int jac_nnz_max;

  int* jac_rows;
  int* jac_cols;
  double* jac_vals;
  int* jac_indices;

  logical goth;

} CUTestConsFuncData;

typedef struct JacCmpData
{
  int* jac_rows;
  int* jac_cols;
} JacCmpData;

static _Thread_local JacCmpData jac_cmp_data;

static int
jac_compare(const void* f, const void* s)
{
  const int first_col  = jac_cmp_data.jac_cols[*((int*)f)];
  const int second_col = jac_cmp_data.jac_cols[*((int*)s)];

  if (first_col < second_col)
  {
    return -1;
  }
  else if (first_col > second_col)
  {
    return 1;
  }

  const int first_row  = jac_cmp_data.jac_rows[*((int*)f)];
  const int second_row = jac_cmp_data.jac_rows[*((int*)s)];

  if (first_row < second_row)
  {
    return -1;
  }
  else if (first_row > second_row)
  {
    return 1;
  }

  return 0;
}

static SLEQP_RETCODE
cutest_cons_data_create(CUTestConsFuncData** star,
                        int num_variables,
                        int num_constraints,
                        int num_linear,
                        SleqpSettings* settings)
{
  SLEQP_CALL(sleqp_malloc(star));

  CUTestConsFuncData* data = *star;
  int status;

  data->eps      = sleqp_settings_real_value(settings, SLEQP_SETTINGS_REAL_EPS);
  data->zero_eps = sleqp_settings_real_value(settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  data->num_constraints = num_constraints;
  data->num_variables   = num_variables;
  data->num_linear      = num_linear;
  data->goth            = cutest_false;

  SLEQP_CALL(sleqp_alloc_array(&data->x, num_variables));
  SLEQP_CALL(sleqp_alloc_array(&data->cons_vals, num_constraints));

  SLEQP_CALL(sleqp_alloc_array(&data->obj_grad, num_variables));

  SLEQP_CALL(sleqp_alloc_array(&data->direction, num_variables));
  SLEQP_CALL(sleqp_alloc_array(&data->multipliers, data->num_constraints));

  SLEQP_CALL(sleqp_alloc_array(&data->hessian_product, num_variables));
  SLEQP_CALL(sleqp_alloc_array(&data->cons_hessian_product, num_variables));

  CUTEST_cdimsj(&status, &(data->jac_nnz_max));

  SLEQP_CUTEST_CHECK_STATUS(status);

  SLEQP_CALL(sleqp_alloc_array(&data->jac_rows, data->jac_nnz_max));
  SLEQP_CALL(sleqp_alloc_array(&data->jac_cols, data->jac_nnz_max));
  SLEQP_CALL(sleqp_alloc_array(&data->jac_vals, data->jac_nnz_max));
  SLEQP_CALL(sleqp_alloc_array(&data->jac_indices, data->jac_nnz_max));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
cutest_cons_data_free(void* func_data)
{
  CUTestConsFuncData* data  = (CUTestConsFuncData*)func_data;
  CUTestConsFuncData** star = &data;

  sleqp_free(&data->jac_indices);
  sleqp_free(&data->jac_vals);
  sleqp_free(&data->jac_cols);
  sleqp_free(&data->jac_rows);

  sleqp_free(&data->cons_hessian_product);
  sleqp_free(&data->hessian_product);

  sleqp_free(&data->multipliers);
  sleqp_free(&data->direction);

  sleqp_free(&data->obj_grad);

  sleqp_free(&data->cons_vals);
  sleqp_free(&data->x);

  sleqp_free(star);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
cutest_cons_func_set(SleqpFunc* func,
                     SleqpVec* x,
                     SLEQP_VALUE_REASON reason,
                     bool* reject,
                     void* func_data)
{
  CUTestConsFuncData* data = (CUTestConsFuncData*)func_data;

  SLEQP_CALL(sleqp_vec_to_raw(x, data->x));

  data->goth = cutest_false;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
cutest_cons_func_nonzeros(SleqpFunc* func,
                          int* obj_grad_nnz,
                          int* cons_val_nnz,
                          int* cons_jac_nnz,
                          int* hess_prod_nnz,
                          void* func_data)
{
  CUTestConsFuncData* data = (CUTestConsFuncData*)func_data;

  *cons_jac_nnz = data->jac_nnz_max;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
cutest_cons_func_set_raw(SleqpFunc* func, const double* x)
{
  void* func_data          = sleqp_func_get_data(func);
  CUTestConsFuncData* data = (CUTestConsFuncData*)func_data;

  const int num_variables = data->num_variables;

  for (int j = 0; j < num_variables; ++j)
  {
    data->x[j] = x[j];
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
cutest_cons_func_obj_val(SleqpFunc* func, double* obj_val, void* func_data)
{
  CUTestConsFuncData* data = (CUTestConsFuncData*)func_data;
  int status;

  CUTEST_cfn(&status,
             &data->num_variables,
             &data->num_constraints,
             data->x,
             obj_val,
             data->cons_vals);

  SLEQP_CUTEST_CHECK_STATUS(status);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
cutest_cons_func_obj_grad(SleqpFunc* func, SleqpVec* obj_grad, void* func_data)
{
  CUTestConsFuncData* data = (CUTestConsFuncData*)func_data;
  int status;

  CUTEST_csgr(&status,                // status flag
              &data->num_variables,   // number of variables
              &data->num_constraints, // number of constraints
              data->x,                // current iterate
              NULL,                   // Lagrangian multipliers
              &cutest_false,    // Do we want the gradient of the Lagrangian?
              &(data->jac_nnz), // Actual number of Jacobian nonzeroes
              &(data->jac_nnz_max), // Maximum number of Jacobian nonzeroes
              data->jac_vals,       // Jacobian data
              data->jac_cols,       // Lagrangian leading size
              data->jac_rows);      // Lagrangian trailing size

  SLEQP_CUTEST_CHECK_STATUS(status);

  for (int j = 0; j < data->num_variables; ++j)
  {
    data->obj_grad[j] = 0.;
  }

  for (int i = 0; i < data->jac_nnz; ++i)
  {
    int row    = data->jac_rows[i];
    int col    = data->jac_cols[i];
    double val = data->jac_vals[i];

    --col;

    assert(col >= 0);

    if (row != 0)
    {
      continue;
    }

    data->obj_grad[col] = val;
  }

  SLEQP_CALL(sleqp_vec_set_from_raw(obj_grad,
                                    data->obj_grad,
                                    data->num_variables,
                                    data->zero_eps));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
cutest_cons_cons_val(SleqpFunc* func, SleqpVec* cons_val, void* func_data)
{
  CUTestConsFuncData* data = (CUTestConsFuncData*)func_data;
  int status;

  double obj;

  const int num_general = data->num_constraints - data->num_linear;

  CUTEST_cfn(&status,
             &data->num_variables,
             &data->num_constraints,
             data->x,
             &obj,
             data->cons_vals);

  SLEQP_CUTEST_CHECK_STATUS(status);

  SLEQP_CALL(sleqp_vec_set_from_raw(cons_val,
                                    data->cons_vals,
                                    num_general,
                                    data->zero_eps));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
cutest_cons_cons_jac(SleqpFunc* func, SleqpMat* cons_jac, void* func_data)
{
  CUTestConsFuncData* data = (CUTestConsFuncData*)func_data;
  int status;

  CUTEST_csgr(&status,                // status flag
              &data->num_variables,   // number of variables
              &data->num_constraints, // number of constraints
              data->x,                // current iterate
              NULL,                   // Lagrangian multipliers
              &cutest_false,    // Do we want the gradient of the Lagrangian?
              &(data->jac_nnz), // Actual number of Jacobian nonzeroes
              &(data->jac_nnz_max), // Maximum number of Jacobian nonzeroes
              data->jac_vals,       // Jacobian data
              data->jac_cols,       // Lagrangian leading size
              data->jac_rows);      // Lagrangian trailing size

  SLEQP_CUTEST_CHECK_STATUS(status);

  for (int i = 0; i < data->jac_nnz; ++i)
  {
    data->jac_indices[i] = i;
  }

  jac_cmp_data.jac_cols = data->jac_cols;
  jac_cmp_data.jac_rows = data->jac_rows;

  qsort(data->jac_indices, data->jac_nnz, sizeof(int), &jac_compare);

  const int num_cols = sleqp_mat_num_cols(cons_jac);
  int last_col       = 0;

  SLEQP_CALL(sleqp_mat_reserve(cons_jac, data->jac_nnz));

  const int num_general = data->num_constraints - data->num_linear;

  for (int i = 0; i < data->jac_nnz; ++i)
  {
    const int k = data->jac_indices[i];

    int row    = data->jac_rows[k];
    int col    = data->jac_cols[k];
    double val = data->jac_vals[k];

    --col;

    assert(col >= 0);

    if (row == 0)
    {
      continue;
    }

    --row;

    // linear constraint
    if (row >= num_general)
    {
      continue;
    }

    while (col > last_col)
    {
      SLEQP_CALL(sleqp_mat_push_col(cons_jac, ++last_col));
    }

    last_col = col;

    SLEQP_CALL(sleqp_mat_push(cons_jac, row, col, val));
  }

  ++last_col;
  while (num_cols > last_col)
  {
    SLEQP_CALL(sleqp_mat_push_col(cons_jac, last_col++));
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_cutest_eval_linear(SleqpFunc* func, SleqpVec* linear)
{
  void* func_data          = sleqp_func_get_data(func);
  CUTestConsFuncData* data = (CUTestConsFuncData*)func_data;
  int status;

  double obj;

  assert(linear->dim == data->num_linear);

  const int num_general = data->num_constraints - data->num_linear;

  CUTEST_cfn(&status,
             &data->num_variables,
             &data->num_constraints,
             data->x,
             &obj,
             data->cons_vals);

  SLEQP_CUTEST_CHECK_STATUS(status);

  SLEQP_CALL(sleqp_vec_set_from_raw(linear,
                                    data->cons_vals + num_general,
                                    data->num_linear,
                                    data->zero_eps));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_cutest_eval_linear_coeffs(SleqpFunc* func, SleqpMat* coeffs)
{
  void* func_data          = sleqp_func_get_data(func);
  CUTestConsFuncData* data = (CUTestConsFuncData*)func_data;
  int status;

  assert(sleqp_mat_num_cols(coeffs) == data->num_variables);
  assert(sleqp_mat_num_rows(coeffs) == data->num_linear);

  CUTEST_csgr(&status,                // status flag
              &data->num_variables,   // number of variables
              &data->num_constraints, // number of constraints
              data->x,                // current iterate
              NULL,                   // Lagrangian multipliers
              &cutest_false,    // Do we want the gradient of the Lagrangian?
              &(data->jac_nnz), // Actual number of Jacobian nonzeroes
              &(data->jac_nnz_max), // Maximum number of Jacobian nonzeroes
              data->jac_vals,       // Jacobian data
              data->jac_cols,       // Lagrangian leading size
              data->jac_rows);      // Lagrangian trailing size

  SLEQP_CUTEST_CHECK_STATUS(status);

  for (int i = 0; i < data->jac_nnz; ++i)
  {
    data->jac_indices[i] = i;
  }

  jac_cmp_data.jac_cols = data->jac_cols;
  jac_cmp_data.jac_rows = data->jac_rows;

  qsort(data->jac_indices, data->jac_nnz, sizeof(int), &jac_compare);

  const int num_cols = sleqp_mat_num_cols(coeffs);
  int last_col       = 0;

  SLEQP_CALL(sleqp_mat_reserve(coeffs, data->jac_nnz));

  const int num_general = data->num_constraints - data->num_linear;

  for (int i = 0; i < data->jac_nnz; ++i)
  {
    const int k = data->jac_indices[i];

    int row    = data->jac_rows[k];
    int col    = data->jac_cols[k];
    double val = data->jac_vals[k];

    --col;

    assert(col >= 0);

    // objective gradient
    if (row == 0)
    {
      continue;
    }

    --row;

    // general constraints
    if (row < num_general)
    {
      continue;
    }

    row -= num_general;

    while (col > last_col)
    {
      SLEQP_CALL(sleqp_mat_push_col(coeffs, ++last_col));
    }

    last_col = col;

    SLEQP_CALL(sleqp_mat_push(coeffs, row, col, val));
  }

  ++last_col;

  while (num_cols > last_col)
  {
    SLEQP_CALL(sleqp_mat_push_col(coeffs, last_col++));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
cutest_cons_func_hess_product(SleqpFunc* func,
                              const SleqpVec* direction,
                              const SleqpVec* cons_duals,
                              SleqpVec* product,
                              void* func_data)
{
  CUTestConsFuncData* data = (CUTestConsFuncData*)func_data;
  int status;

  SLEQP_CALL(sleqp_vec_to_raw(direction, data->direction));

  SLEQP_CALL(sleqp_vec_to_raw(cons_duals, data->multipliers));

  /*
  {
    CUTEST_uhprod(&status,
                  &(data->num_variables),
                  &(data->goth),
                  data->x,
                  data->direction,
                  data->hessian_product);

    SLEQP_CUTEST_CHECK_STATUS(status);

    data->goth = cutest_true;

    CUTEST_chcprod(&status,
                   &(data->num_variables),
                   &(data->num_constraints),
                   &(data->goth),
                   data->x,
                   data->multipliers,
                   data->direction,
                   data->cons_hessian_product);

    SLEQP_CUTEST_CHECK_STATUS(status);

    for(int i = 0; i < data->num_variables; ++i)
    {
      data->hessian_product[i] += data->cons_hessian_product[i];
    }
  }
  */

  {
    CUTEST_chprod(&status,
                  &(data->num_variables),
                  &(data->num_constraints),
                  &(data->goth),
                  data->x,
                  data->multipliers,
                  data->direction,
                  data->hessian_product);

    SLEQP_CUTEST_CHECK_STATUS(status);

    data->goth = cutest_true;
  }

  SLEQP_CALL(sleqp_vec_set_from_raw(product,
                                    data->hessian_product,
                                    data->num_variables,
                                    data->zero_eps));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_cutest_cons_func_create(SleqpFunc** star,
                              int num_variables,
                              int num_constraints,
                              int num_linear,
                              SleqpSettings* settings)
{
  CUTestConsFuncData* data;

  const int num_general = num_constraints - num_linear;

  SLEQP_CALL(cutest_cons_data_create(&data,
                                     num_variables,
                                     num_constraints,
                                     num_linear,
                                     settings));

  SleqpFuncCallbacks callbacks = {.set_value = cutest_cons_func_set,
                                  .nonzeros  = cutest_cons_func_nonzeros,
                                  .obj_val   = cutest_cons_func_obj_val,
                                  .obj_grad  = cutest_cons_func_obj_grad,
                                  .cons_val  = cutest_cons_cons_val,
                                  .cons_jac  = cutest_cons_cons_jac,
                                  .hess_prod = cutest_cons_func_hess_product,
                                  .func_free = cutest_cons_data_free};

  SLEQP_CALL(
    sleqp_func_create(star, &callbacks, num_variables, num_general, data));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
adjust_for_linear_offset(SleqpSettings* settings,
                         SleqpCutestData* data,
                         const SleqpVec* linear_offset,
                         SleqpVec* sparse_cache,
                         SleqpVec* linear_lb,
                         SleqpVec* linear_ub)
{
  const int num_constraints = data->num_constraints;
  const int num_linear      = data->num_linear;
  const int num_general     = num_constraints - num_linear;

  const double zero_eps = sleqp_settings_real_value(settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  SLEQP_CALL(sleqp_vec_set_from_raw(sparse_cache,
                                    data->cons_lb + num_general,
                                    num_linear,
                                    zero_eps));

  SLEQP_CALL(sleqp_vec_add_scaled(sparse_cache,
                                  linear_offset,
                                  1.,
                                  -1.,
                                  zero_eps,
                                  linear_lb));

  SLEQP_CALL(sleqp_vec_set_from_raw(sparse_cache,
                                    data->cons_ub + num_general,
                                    num_linear,
                                    zero_eps));

  SLEQP_CALL(sleqp_vec_add_scaled(sparse_cache,
                                  linear_offset,
                                  1.,
                                  -1.,
                                  zero_eps,
                                  linear_ub));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_linear_offset(SleqpFunc* func,
                      SleqpSettings* settings,
                      SleqpCutestData* data,
                      const SleqpMat* linear_coeffs,
                      SleqpVec* linear_offset)
{
  SleqpVec* linear;
  SleqpVec* x;
  SleqpVec* product;
  double* raw_product;

  const int num_variables = data->num_variables;
  const int num_linear    = data->num_linear;

  const double zero_eps = sleqp_settings_real_value(settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  SLEQP_CALL(sleqp_vec_create_empty(&linear, num_linear));
  SLEQP_CALL(sleqp_vec_create_empty(&x, num_variables));
  SLEQP_CALL(sleqp_vec_create_empty(&product, num_linear));
  SLEQP_CALL(sleqp_alloc_array(&raw_product, num_linear));

  SLEQP_CALL(sleqp_vec_set_from_raw(x, data->x, num_variables, zero_eps));

  SLEQP_CALL(sleqp_cutest_eval_linear(func, linear));

  SLEQP_CALL(sleqp_mat_mult_vec(linear_coeffs, x, raw_product));

  SLEQP_CALL(
    sleqp_vec_set_from_raw(product, raw_product, num_linear, zero_eps));

  SLEQP_CALL(
    sleqp_vec_add_scaled(linear, product, 1., -1., zero_eps, linear_offset));

  sleqp_free(&raw_product);

  SLEQP_CALL(sleqp_vec_free(&product));
  SLEQP_CALL(sleqp_vec_free(&x));
  SLEQP_CALL(sleqp_vec_free(&linear));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_cutest_cons_problem_create(SleqpProblem** star,
                                 SleqpCutestData* data,
                                 SleqpSettings* settings,
                                 bool force_nonlinear)
{
  const int num_variables   = data->num_variables;
  const int num_constraints = data->num_constraints;
  const int num_linear      = force_nonlinear ? 0 : data->num_linear;
  const int num_general     = num_constraints - num_linear;

  const double zero_eps = sleqp_settings_real_value(settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  SleqpVec* var_lb;
  SleqpVec* var_ub;

  SleqpVec* cons_lb;
  SleqpVec* cons_ub;
  SleqpFunc* func;

  SleqpMat* linear_coeffs;
  SleqpVec* sparse_cache;
  SleqpVec* linear_lb;
  SleqpVec* linear_ub;
  SleqpVec* linear_offset;

  SLEQP_CALL(sleqp_vec_create_empty(&var_lb, num_variables));
  SLEQP_CALL(sleqp_vec_create_empty(&var_ub, num_variables));

  SLEQP_CALL(sleqp_vec_create_empty(&cons_lb, num_general));
  SLEQP_CALL(sleqp_vec_create_empty(&cons_ub, num_general));

  SLEQP_CALL(sleqp_vec_create_empty(&sparse_cache, num_linear));

  SLEQP_CALL(sleqp_vec_create_empty(&linear_lb, num_linear));
  SLEQP_CALL(sleqp_vec_create_empty(&linear_ub, num_linear));

  SLEQP_CALL(sleqp_vec_create_empty(&linear_offset, num_linear));

  SLEQP_CALL(
    sleqp_vec_set_from_raw(var_lb, data->var_lb, num_variables, zero_eps));
  SLEQP_CALL(
    sleqp_vec_set_from_raw(var_ub, data->var_ub, num_variables, zero_eps));

  SLEQP_CALL(
    sleqp_vec_set_from_raw(cons_lb, data->cons_lb, num_general, zero_eps));
  SLEQP_CALL(
    sleqp_vec_set_from_raw(cons_ub, data->cons_ub, num_general, zero_eps));

  SLEQP_CALL(sleqp_mat_create(&linear_coeffs, num_linear, num_variables, 0));

  SLEQP_CALL(sleqp_cutest_cons_func_create(&func,
                                           num_variables,
                                           num_constraints,
                                           num_linear,
                                           settings));

  if (num_linear != 0)
  {
    SLEQP_CALL(cutest_cons_func_set_raw(func, data->x));

    SLEQP_CALL(sleqp_cutest_eval_linear_coeffs(func, linear_coeffs));

    SLEQP_CALL(
      compute_linear_offset(func, settings, data, linear_coeffs, linear_offset));

    SLEQP_CALL(adjust_for_linear_offset(settings,
                                        data,
                                        linear_offset,
                                        sparse_cache,
                                        linear_lb,
                                        linear_ub));

    SLEQP_CALL(sleqp_problem_create(star,
                                    func,
                                    settings,
                                    var_lb,
                                    var_ub,
                                    cons_lb,
                                    cons_ub,
                                    linear_coeffs,
                                    linear_lb,
                                    linear_ub));
  }
  else
  {
    SLEQP_CALL(sleqp_problem_create_simple(star,
                                           func,
                                           settings,
                                           var_lb,
                                           var_ub,
                                           cons_lb,
                                           cons_ub));
  }

  SLEQP_CALL(sleqp_func_release(&func));

  SLEQP_CALL(sleqp_vec_free(&linear_offset));

  SLEQP_CALL(sleqp_vec_free(&linear_ub));
  SLEQP_CALL(sleqp_vec_free(&linear_lb));

  SLEQP_CALL(sleqp_mat_release(&linear_coeffs));

  SLEQP_CALL(sleqp_vec_free(&sparse_cache));

  SLEQP_CALL(sleqp_vec_free(&cons_ub));
  SLEQP_CALL(sleqp_vec_free(&cons_lb));

  SLEQP_CALL(sleqp_vec_free(&var_ub));
  SLEQP_CALL(sleqp_vec_free(&var_lb));

  return SLEQP_OKAY;
}
