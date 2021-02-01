#include "quadfunc_fixture.h"

const int num_variables = 2;
const int num_constraints = 0;

SleqpFunc* quadfunc;

SleqpSparseVec* quadfunc_var_lb;
SleqpSparseVec* quadfunc_var_ub;
SleqpSparseVec* quadfunc_cons_lb;
SleqpSparseVec* quadfunc_cons_ub;
SleqpSparseVec* quadfunc_x;

typedef struct SquareFuncData
{
  double* x;
} SquareFuncData;

SquareFuncData* func_data;

static inline double square(double v)
{
  return v*v;
}

SLEQP_RETCODE quadfunc_set(SleqpFunc* func,
                           SleqpSparseVec* x,
                           SLEQP_VALUE_REASON reason,
                           int* func_grad_nnz,
                           int* cons_val_nnz,
                           int* cons_jac_nnz,
                           void* func_data)
{
  *func_grad_nnz = 2;
  *cons_val_nnz = 0;
  *cons_jac_nnz = 0;

  SquareFuncData* data = (SquareFuncData*) func_data;

  data->x[0] = 0;
  data->x[1] = 0;

  int k_x = 0;

  while(k_x < x->nnz)
  {
    data->x[x->indices[k_x]] = x->data[k_x];

    ++k_x;
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE quadfunc_val(SleqpFunc* func,
                            double* func_val,
                            void* func_data)
{
  SquareFuncData* data = (SquareFuncData*) func_data;

  *func_val = square(data->x[0]) + square(data->x[1]);

  return SLEQP_OKAY;
}

SLEQP_RETCODE quadfunc_grad(SleqpFunc* func,
                            SleqpSparseVec* func_grad,
                            void* func_data)
{
  SquareFuncData* data = (SquareFuncData*) func_data;

  assert(func_grad->dim == 2);
  assert(func_grad->nnz_max >= 2);

  func_grad->nnz = 0;

  SLEQP_CALL(sleqp_sparse_vector_push(func_grad,
                                      0,
                                      2.*data->x[0]));

  SLEQP_CALL(sleqp_sparse_vector_push(func_grad,
                                      1,
                                      2.*data->x[1]));

  return SLEQP_OKAY;
}

SLEQP_RETCODE quadfunc_hess_prod(SleqpFunc* func,
                                 const double* func_dual,
                                 const SleqpSparseVec* direction,
                                 const SleqpSparseVec* cons_duals,
                                 SleqpSparseVec* result,
                                 void* func_data)
{
  if(func_dual)
  {
    double total_value = 2.* (*func_dual);

    SLEQP_CALL(sleqp_sparse_vector_copy(direction, result));

    SLEQP_CALL(sleqp_sparse_vector_scale(result, total_value));
  }

  return SLEQP_OKAY;
}

void quadfunc_setup()
{
  ASSERT_CALL(sleqp_malloc(&func_data));

  ASSERT_CALL(sleqp_calloc(&func_data->x, 2));

  SleqpFuncCallbacks callbacks = {
    .set_value = quadfunc_set,
    .func_val  = quadfunc_val,
    .func_grad = quadfunc_grad,
    .cons_val  = NULL,
    .cons_jac  = NULL,
    .hess_prod = quadfunc_hess_prod,
    .func_free = NULL
  };

  ASSERT_CALL(sleqp_func_create(&quadfunc,
                                &callbacks,
                                num_variables,
                                num_constraints,
                                func_data));

  ASSERT_CALL(sleqp_sparse_vector_create_full(&quadfunc_var_lb,
                                              num_variables));

  ASSERT_CALL(sleqp_sparse_vector_push(quadfunc_var_lb, 0, 1.));
  ASSERT_CALL(sleqp_sparse_vector_push(quadfunc_var_lb, 1, 2.));

  ASSERT_CALL(sleqp_sparse_vector_create_full(&quadfunc_var_ub,
                                              num_variables));

  ASSERT_CALL(sleqp_sparse_vector_push(quadfunc_var_ub, 0, 2.));
  ASSERT_CALL(sleqp_sparse_vector_push(quadfunc_var_ub, 1, 3.));

  ASSERT_CALL(sleqp_sparse_vector_create_empty(&quadfunc_cons_lb,
                                               num_constraints));

  ASSERT_CALL(sleqp_sparse_vector_create_empty(&quadfunc_cons_ub,
                                               num_constraints));

  ASSERT_CALL(sleqp_sparse_vector_create_full(&quadfunc_x,
                                              num_variables));


  ASSERT_CALL(sleqp_sparse_vector_push(quadfunc_x, 0, 1.));

  ASSERT_CALL(sleqp_sparse_vector_push(quadfunc_x, 1, 2.));
}

void quadfunc_teardown()
{
  ASSERT_CALL(sleqp_sparse_vector_free(&quadfunc_x));

  ASSERT_CALL(sleqp_sparse_vector_free(&quadfunc_cons_ub));

  ASSERT_CALL(sleqp_sparse_vector_free(&quadfunc_cons_lb));

  ASSERT_CALL(sleqp_sparse_vector_free(&quadfunc_var_ub));

  ASSERT_CALL(sleqp_sparse_vector_free(&quadfunc_var_lb));


  ASSERT_CALL(sleqp_func_release(&quadfunc));

  sleqp_free(&func_data->x);

  sleqp_free(&func_data);
}
