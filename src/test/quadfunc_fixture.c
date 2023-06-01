#include "quadfunc_fixture.h"

const int num_variables   = 2;
const int num_constraints = 0;

SleqpFunc* quadfunc;

SleqpVec* quadfunc_var_lb;
SleqpVec* quadfunc_var_ub;
SleqpVec* quadfunc_cons_lb;
SleqpVec* quadfunc_cons_ub;
SleqpVec* quadfunc_x;

typedef struct SquareFuncData
{
  double* x;
} SquareFuncData;

SquareFuncData* func_data;

static inline double
square(double v)
{
  return v * v;
}

SLEQP_RETCODE
quadfunc_set(SleqpFunc* func,
             SleqpVec* x,
             SLEQP_VALUE_REASON reason,
             bool* reject,
             void* func_data)
{
  SquareFuncData* data = (SquareFuncData*)func_data;

  data->x[0] = 0;
  data->x[1] = 0;

  int k_x = 0;

  while (k_x < x->nnz)
  {
    data->x[x->indices[k_x]] = x->data[k_x];

    ++k_x;
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
quadfunc_obj_val(SleqpFunc* func, double* obj_val, void* func_data)
{
  SquareFuncData* data = (SquareFuncData*)func_data;

  *obj_val = square(data->x[0]) + square(data->x[1]);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
quadfunc_obj_grad(SleqpFunc* func, SleqpVec* obj_grad, void* func_data)
{
  SquareFuncData* data = (SquareFuncData*)func_data;

  assert(obj_grad->dim == 2);
  assert(obj_grad->nnz_max >= 2);

  obj_grad->nnz = 0;

  SLEQP_CALL(sleqp_vec_push(obj_grad, 0, 2. * data->x[0]));

  SLEQP_CALL(sleqp_vec_push(obj_grad, 1, 2. * data->x[1]));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
quadfunc_hess_prod(SleqpFunc* func,
                   const SleqpVec* direction,
                   const SleqpVec* cons_duals,
                   SleqpVec* result,
                   void* func_data)
{
  SLEQP_CALL(sleqp_vec_copy(direction, result));

  SLEQP_CALL(sleqp_vec_scale(result, 2.));

  return SLEQP_OKAY;
}

void
quadfunc_setup()
{
  ASSERT_CALL(sleqp_malloc(&func_data));

  ASSERT_CALL(sleqp_alloc_array(&func_data->x, 2));

  SleqpFuncCallbacks callbacks = {.set_value = quadfunc_set,
                                  .obj_val   = quadfunc_obj_val,
                                  .obj_grad  = quadfunc_obj_grad,
                                  .cons_val  = NULL,
                                  .cons_jac  = NULL,
                                  .hess_prod = quadfunc_hess_prod,
                                  .func_free = NULL};

  ASSERT_CALL(sleqp_func_create(&quadfunc,
                                &callbacks,
                                num_variables,
                                num_constraints,
                                func_data));

  ASSERT_CALL(sleqp_vec_create_full(&quadfunc_var_lb, num_variables));

  ASSERT_CALL(sleqp_vec_push(quadfunc_var_lb, 0, 1.));
  ASSERT_CALL(sleqp_vec_push(quadfunc_var_lb, 1, 2.));

  ASSERT_CALL(sleqp_vec_create_full(&quadfunc_var_ub, num_variables));

  ASSERT_CALL(sleqp_vec_push(quadfunc_var_ub, 0, 2.));
  ASSERT_CALL(sleqp_vec_push(quadfunc_var_ub, 1, 3.));

  ASSERT_CALL(sleqp_vec_create_empty(&quadfunc_cons_lb, num_constraints));

  ASSERT_CALL(sleqp_vec_create_empty(&quadfunc_cons_ub, num_constraints));

  ASSERT_CALL(sleqp_vec_create_full(&quadfunc_x, num_variables));

  ASSERT_CALL(sleqp_vec_push(quadfunc_x, 0, 1.));

  ASSERT_CALL(sleqp_vec_push(quadfunc_x, 1, 2.));
}

void
quadfunc_teardown()
{
  ASSERT_CALL(sleqp_vec_free(&quadfunc_x));

  ASSERT_CALL(sleqp_vec_free(&quadfunc_cons_ub));

  ASSERT_CALL(sleqp_vec_free(&quadfunc_cons_lb));

  ASSERT_CALL(sleqp_vec_free(&quadfunc_var_ub));

  ASSERT_CALL(sleqp_vec_free(&quadfunc_var_lb));

  ASSERT_CALL(sleqp_func_release(&quadfunc));

  sleqp_free(&func_data->x);

  sleqp_free(&func_data);
}
