#include "sleqp_cutest_constrained.h"

#include <stdlib.h>

#include "sleqp_cutest_defs.h"

typedef struct CUTestConsFuncData
{
  double eps;

  int num_constraints;

  double* x;
  double* cons_vals;
  double* func_grad;

  double* direction;
  double* multipliers;
  double* hessian_product;

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
                                                   double eps)
{
  SLEQP_CALL(sleqp_malloc(star));

  CUTestConsFuncData* data = *star;
  int status;

  data->eps = eps;
  data->num_constraints = num_constraints;
  data->goth = cutest_false;

  SLEQP_CALL(sleqp_calloc(&data->x, num_variables));
  SLEQP_CALL(sleqp_calloc(&data->cons_vals, num_constraints));

  SLEQP_CALL(sleqp_calloc(&data->func_grad, num_variables));

  SLEQP_CALL(sleqp_calloc(&data->direction, num_variables));
  SLEQP_CALL(sleqp_calloc(&data->multipliers, data->num_constraints));

  SLEQP_CALL(sleqp_calloc(&data->hessian_product, num_variables));

  CUTEST_cdimsj(&status, &(data->jac_nnz_max));

  SLEQP_CUTEST_CHECK_STATUS(status);

  SLEQP_CALL(sleqp_calloc(&data->jac_rows, data->jac_nnz_max));
  SLEQP_CALL(sleqp_calloc(&data->jac_cols, data->jac_nnz_max));
  SLEQP_CALL(sleqp_calloc(&data->jac_vals, data->jac_nnz_max));
  SLEQP_CALL(sleqp_calloc(&data->jac_indices, data->jac_nnz_max));


  return SLEQP_OKAY;
}

static SLEQP_RETCODE sleqp_cutest_cons_data_free(CUTestConsFuncData** star)
{
  CUTestConsFuncData* data = *star;

  sleqp_free(&data->jac_indices);
  sleqp_free(&data->jac_vals);
  sleqp_free(&data->jac_cols);
  sleqp_free(&data->jac_rows);

  sleqp_free(&data->hessian_product);

  sleqp_free(&data->multipliers);
  sleqp_free(&data->direction);

  sleqp_free(&data->func_grad);

  sleqp_free(&data->cons_vals);
  sleqp_free(&data->x);

  sleqp_free(star);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE sleqp_cutest_cons_func_set(SleqpSparseVec* x,
                                                int num_variables,
                                                int* func_grad_nnz,
                                                int* cons_val_nnz,
                                                int* cons_jac_nnz,
                                                void* func_data)
{
  CUTestConsFuncData* data = (CUTestConsFuncData*) func_data;

  SLEQP_CALL(sleqp_sparse_vector_to_raw(x, data->x));

  data->goth = cutest_false;

  *func_grad_nnz = num_variables;

  *cons_val_nnz = data->num_constraints;

  *cons_jac_nnz = data->jac_nnz_max;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE sleqp_cutest_cons_func_eval(int num_variables,
                                                 SleqpSparseVec* cons_indices,
                                                 double* func_val,
                                                 SleqpSparseVec* func_grad,
                                                 SleqpSparseVec* cons_val,
                                                 SleqpSparseMatrix* cons_jac,
                                                 void* func_data)
{
  CUTestConsFuncData* data = (CUTestConsFuncData*) func_data;
  int status;

  if(func_val || cons_val)
  {
    double obj;


    CUTEST_cfn(&status,
               &num_variables,
               &data->num_constraints,
               data->x,
               &obj,
               data->cons_vals);

    SLEQP_CUTEST_CHECK_STATUS(status);

    if(func_val)
    {
      *func_val = obj;
    }

    if(cons_val)
    {
      SLEQP_CALL(sleqp_sparse_vector_from_raw(cons_val,
                                              data->cons_vals,
                                              data->num_constraints,
                                              data->eps));
    }

  }

  if(func_grad || cons_jac)
  {
    CUTEST_csgr(&status,                // status flag
                &num_variables,         // number of variables
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

    if(func_grad)
    {
      SLEQP_CALL(sleqp_sparse_vector_clear(func_grad));
    }

    int last_col = 0;

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
        if(func_grad)
        {
          SLEQP_CALL(sleqp_sparse_vector_push(func_grad, col, val));
        }

        continue;
      }

      if(!cons_jac)
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

    if(cons_jac)
    {
      ++last_col;
      while(cons_jac->num_cols > last_col)
      {
        SLEQP_CALL(sleqp_sparse_matrix_push_column(cons_jac,
                                                   last_col++));
      }
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE sleqp_cutest_cons_func_hess_product(int num_variables,
                                                         double* func_dual,
                                                         SleqpSparseVec* direction,
                                                         SleqpSparseVec* cons_duals,
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
    CUTEST_chprod(&status,
                  &num_variables,
                  &(data->num_constraints),
                  &(data->goth),
                  data->x,
                  data->multipliers,
                  data->direction,
                  data->hessian_product);

    SLEQP_CUTEST_CHECK_STATUS(status);

  }

  SLEQP_CALL(sleqp_sparse_vector_from_raw(product,
                                          data->hessian_product,
                                          num_variables,
                                          data->eps));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_cutest_cons_func_create(SleqpFunc** star,
                                            int num_variables,
                                            int num_constraints,
                                            double eps)
{
  CUTestConsFuncData* data;

  SLEQP_CALL(sleqp_cutest_cons_data_create(&data,
                                           num_variables,
                                           num_constraints,
                                           eps));

  SLEQP_CALL(sleqp_func_create(star,
                               sleqp_cutest_cons_func_set,
                               sleqp_cutest_cons_func_eval,
                               sleqp_cutest_cons_func_hess_product,
                               num_variables,
                               data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_cutest_cons_func_free(SleqpFunc** star)
{
  SleqpFunc* func = *star;

  CUTestConsFuncData* data = (CUTestConsFuncData*) sleqp_func_get_data(func);

  SLEQP_CALL(sleqp_func_free(star));

  SLEQP_CALL(sleqp_cutest_cons_data_free(&data));

  return SLEQP_OKAY;
}
