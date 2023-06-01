#include "sleqp_cutest_unconstrained.h"

#include "fail.h"
#include "log.h"
#include "mem.h"

#include "sleqp_cutest_types.h"

typedef struct CUTestUnconsFuncData
{
  double eps;

  int num_variables;

  double* x;
  double* obj_grad;

  double* direction;
  double* hessian_product;

  logical goth;

} CUTestUnconsFuncData;

static SLEQP_RETCODE
cutest_uncons_data_create(CUTestUnconsFuncData** star,
                          int num_variables,
                          double eps)
{
  SLEQP_CALL(sleqp_malloc(star));

  CUTestUnconsFuncData* data = *star;

  data->eps           = eps;
  data->num_variables = num_variables;
  data->goth          = cutest_false;

  SLEQP_CALL(sleqp_alloc_array(&data->x, num_variables));

  SLEQP_CALL(sleqp_alloc_array(&data->obj_grad, num_variables));

  SLEQP_CALL(sleqp_alloc_array(&data->direction, num_variables));
  SLEQP_CALL(sleqp_alloc_array(&data->hessian_product, num_variables));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
cutest_uncons_data_free(void* data)
{
  CUTestUnconsFuncData* uncons_data = (CUTestUnconsFuncData*)data;
  CUTestUnconsFuncData** star       = &uncons_data;

  sleqp_free(&uncons_data->hessian_product);
  sleqp_free(&uncons_data->direction);

  sleqp_free(&uncons_data->obj_grad);
  sleqp_free(&uncons_data->x);

  sleqp_free(star);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
cutest_uncons_func_set(SleqpFunc* func,
                       SleqpVec* x,
                       SLEQP_VALUE_REASON reason,
                       bool* reject,
                       void* func_data)
{
  CUTestUnconsFuncData* data = (CUTestUnconsFuncData*)func_data;

  SLEQP_CALL(sleqp_vec_to_raw(x, data->x));

  data->goth = cutest_false;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
cutest_uncons_func_obj_val(SleqpFunc* func, double* obj_val, void* func_data)
{
  CUTestUnconsFuncData* data = (CUTestUnconsFuncData*)func_data;
  int status;

  CUTEST_ufn(&status, &data->num_variables, data->x, obj_val);

  SLEQP_CUTEST_CHECK_STATUS(status);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
cutest_uncons_func_obj_grad(SleqpFunc* func,
                            SleqpVec* obj_grad,
                            void* func_data)
{
  CUTestUnconsFuncData* data = (CUTestUnconsFuncData*)func_data;
  int status;

  CUTEST_ugr(&status, &data->num_variables, data->x, data->obj_grad);

  SLEQP_CUTEST_CHECK_STATUS(status);

  SLEQP_CALL(sleqp_vec_set_from_raw(obj_grad,
                                    data->obj_grad,
                                    data->num_variables,
                                    data->eps));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
cutest_uncons_func_hess_product(SleqpFunc* func,
                                const SleqpVec* direction,
                                const SleqpVec* cons_duals,
                                SleqpVec* product,
                                void* func_data)
{
  CUTestUnconsFuncData* data = (CUTestUnconsFuncData*)func_data;
  int status;

  SLEQP_CALL(sleqp_vec_to_raw(direction, data->direction));

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

  SLEQP_CALL(sleqp_vec_set_from_raw(product,
                                    data->hessian_product,
                                    data->num_variables,
                                    data->eps));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_cutest_uncons_func_create(SleqpFunc** star,
                                int num_variables,
                                SleqpParams* params)
{
  CUTestUnconsFuncData* data;

  const int num_constraints = 0;

  const double zero_eps = sleqp_params_value(params, SLEQP_PARAM_ZERO_EPS);

  SLEQP_CALL(cutest_uncons_data_create(&data, num_variables, zero_eps));

  SleqpFuncCallbacks callbacks = {.set_value = cutest_uncons_func_set,
                                  .obj_val   = cutest_uncons_func_obj_val,
                                  .obj_grad  = cutest_uncons_func_obj_grad,
                                  .cons_val  = NULL,
                                  .cons_jac  = NULL,
                                  .hess_prod = cutest_uncons_func_hess_product,
                                  .func_free = cutest_uncons_data_free};

  SLEQP_CALL(
    sleqp_func_create(star, &callbacks, num_variables, num_constraints, data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_cutest_uncons_problem_create(SleqpProblem** star,
                                   SleqpCutestData* data,
                                   SleqpParams* params)
{
  const int num_variables   = data->num_variables;
  const int num_constraints = data->num_constraints;

  assert(num_constraints == 0);

  const double zero_eps = sleqp_params_value(params, SLEQP_PARAM_ZERO_EPS);

  SleqpVec* var_lb;
  SleqpVec* var_ub;

  SleqpVec* cons_lb;
  SleqpVec* cons_ub;

  SleqpFunc* func;

  SLEQP_CALL(sleqp_vec_create_empty(&var_lb, num_variables));
  SLEQP_CALL(sleqp_vec_create_empty(&var_ub, num_variables));

  SLEQP_CALL(sleqp_vec_create_empty(&cons_lb, num_constraints));
  SLEQP_CALL(sleqp_vec_create_empty(&cons_ub, num_constraints));

  SLEQP_CALL(
    sleqp_vec_set_from_raw(var_lb, data->var_lb, num_variables, zero_eps));
  SLEQP_CALL(
    sleqp_vec_set_from_raw(var_ub, data->var_ub, num_variables, zero_eps));

  SLEQP_CALL(sleqp_cutest_uncons_func_create(&func, num_variables, params));

  SLEQP_CALL(sleqp_problem_create_simple(star,
                                         func,
                                         params,
                                         var_lb,
                                         var_ub,
                                         cons_lb,
                                         cons_ub));

  SLEQP_CALL(sleqp_func_release(&func));

  SLEQP_CALL(sleqp_vec_free(&cons_ub));
  SLEQP_CALL(sleqp_vec_free(&cons_lb));

  SLEQP_CALL(sleqp_vec_free(&var_ub));
  SLEQP_CALL(sleqp_vec_free(&var_lb));

  return SLEQP_OKAY;
}
