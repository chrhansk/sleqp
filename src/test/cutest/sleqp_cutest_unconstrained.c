#include "sleqp_cutest_unconstrained.h"

#include "log.h"
#include "mem.h"

#include "sleqp_cutest_defs.h"

typedef struct CUTestUnconsFuncData
{
  double eps;

  int num_variables;

  double* x;
  double* func_grad;

  double* direction;
  double* hessian_product;

  logical goth;

} CUTestUnconsFuncData;

static SLEQP_RETCODE sleqp_cutest_uncons_data_create(CUTestUnconsFuncData** star,
                                                     int num_variables,
                                                     double eps)
{
  SLEQP_CALL(sleqp_malloc(star));

  CUTestUnconsFuncData* data = *star;

  data->eps = eps;
  data->num_variables = num_variables;
  data->goth = cutest_false;

  SLEQP_CALL(sleqp_alloc_array(&data->x, num_variables));

  SLEQP_CALL(sleqp_alloc_array(&data->func_grad, num_variables));

  SLEQP_CALL(sleqp_alloc_array(&data->direction, num_variables));
  SLEQP_CALL(sleqp_alloc_array(&data->hessian_product, num_variables));

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

static SLEQP_RETCODE sleqp_cutest_uncons_func_set(SleqpFunc* func,
                                                  SleqpSparseVec* x,
                                                  SLEQP_VALUE_REASON reason,
                                                  int* func_grad_nnz,
                                                  int* cons_val_nnz,
                                                  int* cons_jac_nnz,
                                                  void* func_data)
{
  CUTestUnconsFuncData* data = (CUTestUnconsFuncData*) func_data;

  SLEQP_CALL(sleqp_sparse_vector_to_raw(x, data->x));

  data->goth = cutest_false;

  *func_grad_nnz = data->num_variables;

  *cons_val_nnz = 0;
  *cons_jac_nnz = 0;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE sleqp_cutest_uncons_func_val(SleqpFunc* func,
                                                  double* func_val,
                                                  void* func_data)
{
  CUTestUnconsFuncData* data = (CUTestUnconsFuncData*) func_data;
  int status;

  CUTEST_ufn(&status,
             &data->num_variables,
             data->x,
             func_val);

  SLEQP_CUTEST_CHECK_STATUS(status);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE sleqp_cutest_uncons_func_grad(SleqpFunc* func,
                                                   SleqpSparseVec* func_grad,
                                                   void* func_data)
{
  CUTestUnconsFuncData* data = (CUTestUnconsFuncData*) func_data;
  int status;

  CUTEST_ugr(&status,
             &data->num_variables,
             data->x,
             data->func_grad);

  SLEQP_CUTEST_CHECK_STATUS(status);

  SLEQP_CALL(sleqp_sparse_vector_from_raw(func_grad,
                                          data->func_grad,
                                          data->num_variables,
                                          data->eps));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE sleqp_cutest_uncons_func_hess_product(SleqpFunc* func,
                                                           const double* func_dual,
                                                           const SleqpSparseVec* direction,
                                                           const SleqpSparseVec* cons_duals,
                                                           SleqpSparseVec* product,
                                                           void* func_data)
{
  CUTestUnconsFuncData* data = (CUTestUnconsFuncData*) func_data;
  int status;

  SLEQP_CALL(sleqp_sparse_vector_to_raw(direction, data->direction));

  {
    CUTEST_uhprod(&status,
                  &(data->num_variables),
                  &(data->goth),
                  data->x,
                  data->direction,
                  data->hessian_product);

    SLEQP_CUTEST_CHECK_STATUS(status);

    data->goth = cutest_true;
  }

  SLEQP_CALL(sleqp_sparse_vector_from_raw(product,
                                          data->hessian_product,
                                          data->num_variables,
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

  SleqpFuncCallbacks callbacks = {
    .set_value = sleqp_cutest_uncons_func_set,
    .func_val = sleqp_cutest_uncons_func_val,
    .func_grad = sleqp_cutest_uncons_func_grad,
    .cons_val = NULL,
    .cons_jac = NULL,
    .hess_prod = sleqp_cutest_uncons_func_hess_product,
    .func_free = NULL
  };

  SLEQP_CALL(sleqp_func_create(star,
                               &callbacks,
                               num_variables,
                               0,
                               data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_cutest_uncons_func_free(SleqpFunc** star)
{
  SleqpFunc* func = *star;

  CUTestUnconsFuncData* data = (CUTestUnconsFuncData*) sleqp_func_get_data(func);

  SLEQP_CALL(sleqp_func_release(star));

  SLEQP_CALL(sleqp_cutest_uncons_data_free(&data));

  return SLEQP_OKAY;
}
