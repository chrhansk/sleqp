#include "sleqp_cutest_unconstrained.h"

#include <assert.h>

#include "log.h"
#include "mem.h"

#include "sleqp_cutest_types.h"

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

static
SLEQP_RETCODE cutest_uncons_data_create(CUTestUnconsFuncData** star,
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

static
SLEQP_RETCODE cutest_uncons_data_free(void* data)
{
  CUTestUnconsFuncData* uncons_data = (CUTestUnconsFuncData*) data;
  CUTestUnconsFuncData** star = &uncons_data;

  sleqp_free(&uncons_data->hessian_product);
  sleqp_free(&uncons_data->direction);

  sleqp_free(&uncons_data->func_grad);
  sleqp_free(&uncons_data->x);

  sleqp_free(star);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
cutest_uncons_func_set(SleqpFunc* func,
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

static
SLEQP_RETCODE cutest_uncons_func_val(SleqpFunc* func,
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

static
SLEQP_RETCODE cutest_uncons_func_grad(SleqpFunc* func,
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

static
SLEQP_RETCODE cutest_uncons_func_hess_product(SleqpFunc* func,
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
                                              SleqpParams* params)
{
  CUTestUnconsFuncData* data;

  const int num_constraints = 0;

  const double zero_eps = sleqp_params_get(params, SLEQP_PARAM_ZERO_EPS);

  SLEQP_CALL(cutest_uncons_data_create(&data,
                                       num_variables,
                                       zero_eps));

  SleqpFuncCallbacks callbacks = {
    .set_value = cutest_uncons_func_set,
    .func_val  = cutest_uncons_func_val,
    .func_grad = cutest_uncons_func_grad,
    .cons_val  = NULL,
    .cons_jac  = NULL,
    .hess_prod = cutest_uncons_func_hess_product,
    .func_free = cutest_uncons_data_free
  };

  SLEQP_CALL(sleqp_func_create(star,
                               &callbacks,
                               num_variables,
                               num_constraints,
                               data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_cutest_uncons_problem_create(SleqpProblem** star,
                                                 SleqpCutestData* data,
                                                 SleqpParams* params)
{
  const int num_variables = data->num_variables;
  const int num_constraints = data->num_constraints;

  assert(num_constraints == 0);

  const double zero_eps = sleqp_params_get(params, SLEQP_PARAM_ZERO_EPS);

  SleqpSparseVec* var_lb;
  SleqpSparseVec* var_ub;

  SleqpSparseVec* cons_lb;
  SleqpSparseVec* cons_ub;

  SleqpFunc* func;

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&var_lb, num_variables));
  SLEQP_CALL(sleqp_sparse_vector_create_empty(&var_ub, num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&cons_lb, num_constraints));
  SLEQP_CALL(sleqp_sparse_vector_create_empty(&cons_ub, num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_from_raw(var_lb, data->var_lb, num_variables, zero_eps));
  SLEQP_CALL(sleqp_sparse_vector_from_raw(var_ub, data->var_ub, num_variables, zero_eps));

  SLEQP_CALL(sleqp_cutest_uncons_func_create(&func,
                                             num_variables,
                                             params));

  SLEQP_CALL(sleqp_problem_create_simple(star,
                                         func,
                                         params,
                                         var_lb,
                                         var_ub,
                                         cons_lb,
                                         cons_ub));

  SLEQP_CALL(sleqp_func_release(&func));

  SLEQP_CALL(sleqp_sparse_vector_free(&cons_ub));
  SLEQP_CALL(sleqp_sparse_vector_free(&cons_lb));

  SLEQP_CALL(sleqp_sparse_vector_free(&var_ub));
  SLEQP_CALL(sleqp_sparse_vector_free(&var_lb));


  return SLEQP_OKAY;
}
