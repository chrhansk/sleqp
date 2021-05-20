#include "sleqp_cutest_constrained.h"

#include <stdlib.h>

#include "sleqp_cutest_defs.h"

typedef struct CUTestConsFuncData
{
  double eps;
  double zero_eps;

  int num_variables;
  int num_constraints;

  double* x;
  double* cons_vals;
  double* func_grad;

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

static int jac_compare(const void* f, const void* s)
{
  const int first_col = jac_cmp_data.jac_cols[*( (int*) f)];
  const int second_col = jac_cmp_data.jac_cols[*( (int*) s)];

  if(first_col < second_col)
  {
    return -1;
  }
  else if(first_col > second_col)
  {
    return 1;
  }

  const int first_row = jac_cmp_data.jac_rows[*( (int*) f)];
  const int second_row = jac_cmp_data.jac_rows[*( (int*) s)];

  if(first_row < second_row)
  {
    return -1;
  }
  else if(first_row > second_row)
  {
    return 1;
  }

  return 0;
}

static SLEQP_RETCODE sleqp_cutest_cons_data_create(CUTestConsFuncData** star,
                                                   int num_variables,
                                                   int num_constraints,
                                                   SleqpParams* params)
{
  SLEQP_CALL(sleqp_malloc(star));

  CUTestConsFuncData* data = *star;
  int status;

  data->eps = sleqp_params_get(params, SLEQP_PARAM_EPS);
  data->zero_eps = sleqp_params_get(params, SLEQP_PARAM_ZERO_EPS);

  data->num_constraints = num_constraints;
  data->num_variables = num_variables;
  data->goth = cutest_false;

  SLEQP_CALL(sleqp_alloc_array(&data->x, num_variables));
  SLEQP_CALL(sleqp_alloc_array(&data->cons_vals, num_constraints));

  SLEQP_CALL(sleqp_alloc_array(&data->func_grad, num_variables));

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

static SLEQP_RETCODE sleqp_cutest_cons_data_free(CUTestConsFuncData** star)
{
  CUTestConsFuncData* data = *star;

  sleqp_free(&data->jac_indices);
  sleqp_free(&data->jac_vals);
  sleqp_free(&data->jac_cols);
  sleqp_free(&data->jac_rows);

  sleqp_free(&data->cons_hessian_product);
  sleqp_free(&data->hessian_product);

  sleqp_free(&data->multipliers);
  sleqp_free(&data->direction);

  sleqp_free(&data->func_grad);

  sleqp_free(&data->cons_vals);
  sleqp_free(&data->x);

  sleqp_free(star);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE sleqp_cutest_cons_func_set(SleqpFunc* func,
                                                SleqpSparseVec* x,
                                                SLEQP_VALUE_REASON reason,
                                                int* func_grad_nnz,
                                                int* cons_val_nnz,
                                                int* cons_jac_nnz,
                                                void* func_data)
{
  CUTestConsFuncData* data = (CUTestConsFuncData*) func_data;

  SLEQP_CALL(sleqp_sparse_vector_to_raw(x, data->x));

  data->goth = cutest_false;

  *func_grad_nnz = data->num_variables;

  *cons_val_nnz = data->num_constraints;

  *cons_jac_nnz = data->jac_nnz_max;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE sleqp_cutest_cons_func_val(SleqpFunc* func,
                                                double* func_val,
                                                void* func_data)
{
  CUTestConsFuncData* data = (CUTestConsFuncData*) func_data;
  int status;

  CUTEST_cfn(&status,
             &data->num_variables,
             &data->num_constraints,
             data->x,
             func_val,
             data->cons_vals);

  SLEQP_CUTEST_CHECK_STATUS(status);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE sleqp_cutest_cons_func_grad(SleqpFunc* func,
                                                 SleqpSparseVec* func_grad,
                                                 void* func_data)
{
  CUTestConsFuncData* data = (CUTestConsFuncData*) func_data;
  int status;

  CUTEST_ugr(&status,                // status flag
             &data->num_variables,   // number of variables
             data->x,                // current iterate
             data->func_grad);       // function gradient

  SLEQP_CUTEST_CHECK_STATUS(status);

  SLEQP_CALL(sleqp_sparse_vector_from_raw(func_grad,
                                          data->func_grad,
                                          data->num_variables,
                                          data->zero_eps));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE sleqp_cutest_cons_cons_val(SleqpFunc* func,
                                                const SleqpSparseVec* cons_indices,
                                                SleqpSparseVec* cons_val,
                                                void* func_data)
{
  CUTestConsFuncData* data = (CUTestConsFuncData*) func_data;
  int status;

  double obj;

  CUTEST_cfn(&status,
             &data->num_variables,
             &data->num_constraints,
             data->x,
             &obj,
             data->cons_vals);

  SLEQP_CUTEST_CHECK_STATUS(status);

  SLEQP_CALL(sleqp_sparse_vector_from_raw(cons_val,
                                          data->cons_vals,
                                          data->num_constraints,
                                          data->zero_eps));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE sleqp_cutest_cons_cons_jac(SleqpFunc* func,
                                                const SleqpSparseVec* cons_indices,
                                                SleqpSparseMatrix* cons_jac,
                                                void* func_data)
{
  CUTestConsFuncData* data = (CUTestConsFuncData*) func_data;
  int status;

  CUTEST_csgr(&status,                // status flag
              &data->num_variables,   // number of variables
              &data->num_constraints, // number of constraints
              data->x,                // current iterate
              NULL,                   // Lagrangian multipliers
              &cutest_false,          // Do we want the gradient of the Lagrangian?
              &(data->jac_nnz),       // Actual number of Jacobian nonzeroes
              &(data->jac_nnz_max),   // Maximum number of Jacobian nonzeroes
              data->jac_vals,         // Jacobian data
              data->jac_cols,         // Lagrangian leading size
              data->jac_rows);        // Lagrangian trailing size


  SLEQP_CUTEST_CHECK_STATUS(status);

  for(int i = 0; i < data->jac_nnz; ++i)
  {
    data->jac_indices[i] = i;
  }

  //JacCmpData jac_data = { data->jac_rows, data->jac_cols};

  jac_cmp_data.jac_cols = data->jac_cols;
  jac_cmp_data.jac_rows = data->jac_rows;

  qsort(data->jac_indices,
        data->jac_nnz,
        sizeof(int),
        &jac_compare);

  const int num_cols = sleqp_sparse_matrix_get_num_cols(cons_jac);
  int last_col = 0;

  SLEQP_CALL(sleqp_sparse_matrix_reserve(cons_jac, data->jac_nnz));

  for(int i = 0; i < data->jac_nnz; ++i)
  {
    const int k = data->jac_indices[i];

    int row = data->jac_rows[k];
    int col = data->jac_cols[k];
    double val = data->jac_vals[k];

    --col;

    assert(col >= 0);

    if(row == 0)
    {
      continue;
    }

    --row;

    while(col > last_col)
    {
      SLEQP_CALL(sleqp_sparse_matrix_push_column(cons_jac,
                                                 ++last_col));
    }

    last_col = col;

    //sleqp_log_debug("Pushing row = %d, col = %d, val = %f", row, col, val);

    SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac, row, col, val));
  }

  ++last_col;
  while(num_cols > last_col)
  {
    SLEQP_CALL(sleqp_sparse_matrix_push_column(cons_jac,
                                               last_col++));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE sleqp_cutest_cons_func_hess_product(SleqpFunc* func,
                                                         const double* func_dual,
                                                         const SleqpSparseVec* direction,
                                                         const SleqpSparseVec* cons_duals,
                                                         SleqpSparseVec* product,
                                                         void* func_data)
{
  CUTestConsFuncData* data = (CUTestConsFuncData*) func_data;
  int status;

  assert(func_dual);
  assert(*func_dual == 1.);

  SLEQP_CALL(sleqp_sparse_vector_to_raw(direction, data->direction));

  SLEQP_CALL(sleqp_sparse_vector_to_raw(cons_duals, data->multipliers));

  {
    CUTEST_uhprod(&status,
                  &(data->num_variables),
                  &(data->goth),
                  data->x,
                  data->direction,
                  data->hessian_product);

    SLEQP_CUTEST_CHECK_STATUS(status);

    CUTEST_chcprod(&status,
                   &(data->num_variables),
                   &(data->num_constraints),
                   &(data->goth),
                   data->x,
                   data->multipliers,
                   data->direction,
                   data->cons_hessian_product);

    SLEQP_CUTEST_CHECK_STATUS(status);
  }

  for(int i = 0; i < data->num_variables; ++i)
  {
    data->hessian_product[i] += data->cons_hessian_product[i];
  }

  SLEQP_CALL(sleqp_sparse_vector_from_raw(product,
                                          data->hessian_product,
                                          data->num_variables,
                                          data->zero_eps));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_cutest_cons_func_create(SleqpFunc** star,
                                            int num_variables,
                                            int num_constraints,
                                            SleqpParams* params)
{
  CUTestConsFuncData* data;

  SLEQP_CALL(sleqp_cutest_cons_data_create(&data,
                                           num_variables,
                                           num_constraints,
                                           params));

  SleqpFuncCallbacks callbacks = {
    .set_value = sleqp_cutest_cons_func_set,
    .func_val = sleqp_cutest_cons_func_val,
    .func_grad = sleqp_cutest_cons_func_grad,
    .cons_val = sleqp_cutest_cons_cons_val,
    .cons_jac = sleqp_cutest_cons_cons_jac,
    .hess_prod = sleqp_cutest_cons_func_hess_product,
    .func_free = NULL
  };

  SLEQP_CALL(sleqp_func_create(star,
                               &callbacks,
                               num_variables,
                               num_constraints,
                               data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_cutest_cons_func_free(SleqpFunc** star)
{
  SleqpFunc* func = *star;

  CUTestConsFuncData* data = (CUTestConsFuncData*) sleqp_func_get_data(func);

  SLEQP_CALL(sleqp_func_release(star));

  SLEQP_CALL(sleqp_cutest_cons_data_free(&data));

  return SLEQP_OKAY;
}
