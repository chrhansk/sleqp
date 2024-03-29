#include "constrained_fixture.h"

#include "mem.h"

// Note: This problem is cutest problem HS71

int constrained_num_variables   = 4;
int constrained_num_constraints = 2;

SleqpVec* constrained_var_lb  = NULL;
SleqpVec* constrained_var_ub  = NULL;
SleqpVec* constrained_cons_lb = NULL;
SleqpVec* constrained_cons_ub = NULL;
SleqpVec* constrained_initial = NULL;
SleqpVec* constrained_optimum = NULL;

SleqpFunc* constrained_func = NULL;

typedef struct FuncData
{
  double* values;
  double* duals;
  double* direction;

} FuncData;

static double
sq(double x)
{
  return x * x;
}

static SLEQP_RETCODE
func_set(SleqpFunc* func,
         SleqpVec* x,
         SLEQP_VALUE_REASON reason,
         bool* reject,
         void* func_data)
{
  FuncData* data = (FuncData*)func_data;

  SLEQP_CALL(sleqp_vec_to_raw(x, data->values));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
obj_val(SleqpFunc* func, double* obj_val, void* func_data)
{
  FuncData* data = (FuncData*)func_data;
  double* x      = data->values;

  (*obj_val) = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2];

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
obj_grad(SleqpFunc* func, SleqpVec* obj_grad, void* func_data)
{
  FuncData* data = (FuncData*)func_data;
  double* x      = data->values;

  SLEQP_CALL(
    sleqp_vec_push(obj_grad, 0, (x[0] + x[1] + x[2]) * x[3] + x[0] * x[3]));

  SLEQP_CALL(sleqp_vec_push(obj_grad, 1, x[0] * x[3]));

  SLEQP_CALL(sleqp_vec_push(obj_grad, 2, x[0] * x[3] + 1));

  SLEQP_CALL(sleqp_vec_push(obj_grad, 3, (x[0] + x[1] + x[2]) * x[0]));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
cons_val(SleqpFunc* func, SleqpVec* cons_val, void* func_data)
{
  FuncData* data = (FuncData*)func_data;
  double* x      = data->values;

  SLEQP_CALL(sleqp_vec_push(cons_val, 0, x[0] * x[1] * x[2] * x[3]));

  const double sq_sum = sq(x[0]) + sq(x[1]) + sq(x[2]) + sq(x[3]);

  SLEQP_CALL(sleqp_vec_push(cons_val, 1, sq_sum));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
cons_jac(SleqpFunc* func, SleqpMat* cons_jac, void* func_data)
{
  FuncData* data = (FuncData*)func_data;
  double* x      = data->values;

  SLEQP_CALL(sleqp_mat_push(cons_jac, 0, 0, x[1] * x[2] * x[3]));

  SLEQP_CALL(sleqp_mat_push(cons_jac, 1, 0, 2 * x[0]));

  SLEQP_CALL(sleqp_mat_push_col(cons_jac, 1));

  SLEQP_CALL(sleqp_mat_push(cons_jac, 0, 1, x[0] * x[2] * x[3]));

  SLEQP_CALL(sleqp_mat_push(cons_jac, 1, 1, 2 * x[1]));

  SLEQP_CALL(sleqp_mat_push_col(cons_jac, 2));

  SLEQP_CALL(sleqp_mat_push(cons_jac, 0, 2, x[0] * x[1] * x[3]));

  SLEQP_CALL(sleqp_mat_push(cons_jac, 1, 2, 2 * x[2]));

  SLEQP_CALL(sleqp_mat_push_col(cons_jac, 3));

  SLEQP_CALL(sleqp_mat_push(cons_jac, 0, 3, x[0] * x[1] * x[2]));

  SLEQP_CALL(sleqp_mat_push(cons_jac, 1, 3, 2 * x[3]));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
func_hess_prod(SleqpFunc* func,
               const SleqpVec* direction,
               const SleqpVec* cons_duals,
               SleqpVec* product,
               void* func_data)
{
  FuncData* data = (FuncData*)func_data;

  SLEQP_CALL(sleqp_vec_to_raw(cons_duals, data->duals));

  SLEQP_CALL(sleqp_vec_to_raw(direction, data->direction));

  double* x = data->values;

  double* duals = data->duals;
  double* dir   = data->direction;

  SLEQP_CALL(sleqp_vec_reserve(product, constrained_num_variables));

  {
    double v = 0.;

    v += (2 * x[3] + 2 * duals[1]) * dir[0];
    v += (x[3] + x[2] * x[3] * duals[0]) * dir[1];
    v += (x[3] + x[1] * x[3] * duals[0]) * dir[2];
    v += ((2 * x[0] + x[1] + x[2]) + x[1] * x[2] * duals[0]) * dir[3];

    SLEQP_CALL(sleqp_vec_push(product, 0, v));
  }

  {
    double v = 0.;

    v += (x[3] + x[2] * x[3] * duals[0]) * dir[0];
    v += (2 * duals[1]) * dir[1];
    v += (x[0] * x[3] * duals[0]) * dir[2];
    v += (x[0] + (x[0] * x[2]) * duals[0]) * dir[3];

    SLEQP_CALL(sleqp_vec_push(product, 1, v));
  }

  {
    double v = 0.;

    v += (x[3] + x[1] * x[3] * duals[0]) * dir[0];
    v += (x[0] * x[3] * duals[0]) * dir[1];
    v += (2 * duals[1]) * dir[2];
    v += (x[0] + x[0] * x[1] * duals[0]) * dir[3];

    SLEQP_CALL(sleqp_vec_push(product, 2, v));
  }

  {
    double v = 0.;

    v += ((2 * x[0] + x[1] + x[2]) + x[1] * x[2] * duals[0]) * dir[0];
    v += (x[0] + x[0] * x[2] * duals[0]) * dir[1];
    v += (x[0] + x[0] * x[1] * duals[0]) * dir[2];
    v += (2 * duals[1]) * dir[3];

    SLEQP_CALL(sleqp_vec_push(product, 3, v));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
func_free(void* data)
{
  FuncData* func_data = (FuncData*)data;

  sleqp_free(&func_data->direction);
  sleqp_free(&func_data->duals);
  sleqp_free(&func_data->values);
  sleqp_free(&func_data);

  return SLEQP_OKAY;
}

void
constrained_setup()
{
  const double inf = sleqp_infinity();

  ASSERT_CALL(sleqp_vec_create(&constrained_var_lb,
                               constrained_num_variables,
                               constrained_num_variables));

  ASSERT_CALL(sleqp_vec_create(&constrained_var_ub,
                               constrained_num_variables,
                               constrained_num_variables));

  for (int i = 0; i < constrained_num_variables; ++i)
  {
    ASSERT_CALL(sleqp_vec_push(constrained_var_lb, i, 1.));
    ASSERT_CALL(sleqp_vec_push(constrained_var_ub, i, 5.));
  }

  ASSERT_CALL(sleqp_vec_create(&constrained_cons_lb,
                               constrained_num_constraints,
                               constrained_num_constraints));

  ASSERT_CALL(sleqp_vec_create(&constrained_cons_ub,
                               constrained_num_constraints,
                               constrained_num_constraints));

  ASSERT_CALL(sleqp_vec_push(constrained_cons_lb, 0, 25.));
  ASSERT_CALL(sleqp_vec_push(constrained_cons_lb, 1, 40.));

  ASSERT_CALL(sleqp_vec_push(constrained_cons_ub, 0, inf));
  ASSERT_CALL(sleqp_vec_push(constrained_cons_ub, 1, 40.));

  ASSERT_CALL(sleqp_vec_create(&constrained_initial,
                               constrained_num_variables,
                               constrained_num_variables));

  ASSERT_CALL(sleqp_vec_push(constrained_initial, 0, 1.));
  ASSERT_CALL(sleqp_vec_push(constrained_initial, 1, 5.));
  ASSERT_CALL(sleqp_vec_push(constrained_initial, 2, 5.));
  ASSERT_CALL(sleqp_vec_push(constrained_initial, 3, 1.));

  FuncData* func_data;

  ASSERT_CALL(sleqp_malloc(&func_data));

  ASSERT_CALL(sleqp_alloc_array(&func_data->values, constrained_num_variables));
  ASSERT_CALL(
    sleqp_alloc_array(&func_data->duals, constrained_num_constraints));
  ASSERT_CALL(
    sleqp_alloc_array(&func_data->direction, constrained_num_variables));

  SleqpFuncCallbacks callbacks = {.set_value = func_set,
                                  .obj_val   = obj_val,
                                  .obj_grad  = obj_grad,
                                  .cons_val  = cons_val,
                                  .cons_jac  = cons_jac,
                                  .hess_prod = func_hess_prod,
                                  .func_free = func_free};

  ASSERT_CALL(sleqp_func_create(&constrained_func,
                                &callbacks,
                                constrained_num_variables,
                                constrained_num_constraints,
                                func_data));

  ASSERT_CALL(sleqp_vec_create(&constrained_optimum, 4, 4));

  ASSERT_CALL(sleqp_vec_push(constrained_optimum, 0, 1.));
  ASSERT_CALL(sleqp_vec_push(constrained_optimum, 1, 4.742999));
  ASSERT_CALL(sleqp_vec_push(constrained_optimum, 2, 3.821151));
  ASSERT_CALL(sleqp_vec_push(constrained_optimum, 3, 1.379408));
}

void
constrained_teardown()
{
  ASSERT_CALL(sleqp_vec_free(&constrained_optimum));

  ASSERT_CALL(sleqp_func_release(&constrained_func));

  ASSERT_CALL(sleqp_vec_free(&constrained_initial));

  ASSERT_CALL(sleqp_vec_free(&constrained_cons_ub));
  ASSERT_CALL(sleqp_vec_free(&constrained_cons_lb));
  ASSERT_CALL(sleqp_vec_free(&constrained_var_ub));
  ASSERT_CALL(sleqp_vec_free(&constrained_var_lb));
}
