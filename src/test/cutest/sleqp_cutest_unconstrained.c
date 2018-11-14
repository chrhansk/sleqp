#include "sleqp_cutest_unconstrained.h"

#include <cutest.h>

#include "sleqp_cutest_defs.h"

typedef struct CUTestUnconsFuncData
{
  double eps;

  double* x;
  double* func_grad;

  double* direction;
  double* hessian_product;

  bool goth;

} CUTestUnconsFuncData;

static SLEQP_RETCODE sleqp_cutest_uncons_data_create(CUTestUnconsFuncData** star,
                                                     int num_variables,
                                                     double eps)
{
  SLEQP_CALL(sleqp_malloc(star));

  CUTestUnconsFuncData* data = *star;

  data->eps = eps;
  data->goth = false;

  SLEQP_CALL(sleqp_calloc(&data->x, num_variables));

  SLEQP_CALL(sleqp_calloc(&data->func_grad, num_variables));

  SLEQP_CALL(sleqp_calloc(&data->direction, num_variables));
  SLEQP_CALL(sleqp_calloc(&data->hessian_product, num_variables));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE sleqp_cutest_uncons_data_free(CUTestUnconsFuncData** star)
{
  CUTestUnconsFuncData* data = *star;

  sleqp_free(&data->hessian_product);
  sleqp_free(&data->direction);

  sleqp_free(&data->func_grad);
  sleqp_free(&data->x);

  sleqp_free(star);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE sleqp_cutest_uncons_func_set(SleqpSparseVec* x,
                                                  int num_variables,
                                                  int* func_grad_nnz,
                                                  int* cons_val_nnz,
                                                  int* cons_jac_nnz,
                                                  void* func_data)
{
  CUTestUnconsFuncData* data = (CUTestUnconsFuncData*) func_data;

  SLEQP_CALL(sleqp_sparse_vector_to_raw(x, data->x));

  data->goth = false;

  *func_grad_nnz = num_variables;

  *cons_val_nnz = 0;
  *cons_jac_nnz = 0;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE sleqp_cutest_uncons_func_eval(int num_variables,
                                                   SleqpSparseVec* cons_indices,
                                                   double* func_val,
                                                   SleqpSparseVec* func_grad,
                                                   SleqpSparseVec* cons_val,
                                                   SleqpSparseMatrix* cons_jac,
                                                   void* func_data)
{
  CUTestUnconsFuncData* data = (CUTestUnconsFuncData*) func_data;
  int status;

  if(func_val)
  {
    CUTEST_ufn(&status,
               &num_variables,
               data->x,
               func_val);

    SLEQP_CUTEST_CHECK_STATUS(status);
  }

  if(func_grad)
  {
    CUTEST_ugr(&status,
               &num_variables,
               data->x,
               data->func_grad);

    SLEQP_CUTEST_CHECK_STATUS(status);

    SLEQP_CALL(sleqp_sparse_vector_from_raw(func_grad,
                                            data->func_grad,
                                            num_variables,
                                            data->eps));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE sleqp_cutest_uncons_func_hess_product(int num_variables,
                                                           double* func_dual,
                                                           SleqpSparseVec* direction,
                                                           SleqpSparseVec* cons_duals,
                                                           SleqpSparseVec* product,
                                                           void* func_data)
{
  CUTestUnconsFuncData* data = (CUTestUnconsFuncData*) func_data;
  int status;

  SLEQP_CALL(sleqp_sparse_vector_to_raw(direction, data->direction));

  {
    CUTEST_uhprod(&status,
                  &num_variables,
                  &(data->goth),
                  data->x,
                  data->direction,
                  data->hessian_product);

    SLEQP_CUTEST_CHECK_STATUS(status);

    data->goth = true;
  }

  SLEQP_CALL(sleqp_sparse_vector_from_raw(product,
                                          data->hessian_product,
                                          num_variables,
                                          data->eps));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_cutest_uncons_func_create(SleqpFunc** star,
                                              int num_variables,
                                              double eps)
{
  CUTestUnconsFuncData* data;

  SLEQP_CALL(sleqp_cutest_uncons_data_create(&data,
                                             num_variables,
                                             eps));

  SLEQP_CALL(sleqp_func_create(star,
                               sleqp_cutest_uncons_func_set,
                               sleqp_cutest_uncons_func_eval,
                               sleqp_cutest_uncons_func_hess_product,
                               num_variables,
                               data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_cutest_uncons_func_free(SleqpFunc** star)
{
  SleqpFunc* func = *star;

  CUTestUnconsFuncData* data = (CUTestUnconsFuncData*) sleqp_func_get_data(func);

  SLEQP_CALL(sleqp_func_free(star));

  SLEQP_CALL(sleqp_cutest_uncons_data_free(&data));

  return SLEQP_OKAY;
}
