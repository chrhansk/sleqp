#include "quadcons_fixture.h"

SleqpFunc* quadconsfunc;

SleqpSparseVec* quadconsfunc_var_lb;
SleqpSparseVec* quadconsfunc_var_ub;
SleqpSparseVec* quadconsfunc_cons_lb;
SleqpSparseVec* quadconsfunc_cons_ub;
SleqpSparseVec* quadconsfunc_x;

static const int num_variables = 2;
static const int num_constraints = 2;

const int quadconsfunc_num_variables = num_variables;
const int quadconsfunc_num_constraints = num_constraints;

typedef struct SquareFuncData
{
  double* x;
} SquareFuncData;

SquareFuncData* func_data;

static inline double square(double v)
{
  return v*v;
}

SLEQP_RETCODE quadconsfunc_set(SleqpFunc* func,
                               SleqpSparseVec* x,
                               SLEQP_VALUE_REASON reason,
                               bool* reject,
                               int* func_grad_nnz,
                               int* cons_val_nnz,
                               int* cons_jac_nnz,
                               void* func_data)
{
  *func_grad_nnz = 2;
  *cons_val_nnz = 2;
  *cons_jac_nnz = 4;

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

SLEQP_RETCODE quadconsfunc_val(SleqpFunc* func,
                               double* func_val,
                               void* func_data)
{
  SquareFuncData* data = (SquareFuncData*) func_data;

  *func_val = square(data->x[0]) + square(data->x[1]);

  return SLEQP_OKAY;
}

SLEQP_RETCODE quadconsfunc_grad(SleqpFunc* func,
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

SLEQP_RETCODE quadconsfunc_cons_val(SleqpFunc* func,
                                    const SleqpSparseVec* cons_indices,
                                    SleqpSparseVec* cons_val,
                                    void* func_data)
{
  SquareFuncData* data = (SquareFuncData*) func_data;

  assert(cons_val->dim == 2);
  assert(cons_val->nnz_max >= 2);

  cons_val->nnz = 0;

  SLEQP_CALL(sleqp_sparse_vector_push(cons_val,
                                      0,
                                      square(data->x[0]) + square(data->x[1])));

  SLEQP_CALL(sleqp_sparse_vector_push(cons_val,
                                      1,
                                      (square(1 - data->x[0])) + square(1 - data->x[1])));

  return SLEQP_OKAY;
}

SLEQP_RETCODE quadconsfunc_cons_jac(SleqpFunc* func,
                                    const SleqpSparseVec* cons_indices,
                                    SleqpSparseMatrix* cons_jac,
                                    void* func_data)
{
  SquareFuncData* data = (SquareFuncData*) func_data;

  assert(sleqp_sparse_matrix_get_num_rows(cons_jac) == 2);
  assert(sleqp_sparse_matrix_get_num_cols(cons_jac) == 2);
  assert(sleqp_sparse_matrix_get_nnz_max(cons_jac) >= 4);

  SLEQP_CALL(sleqp_sparse_matrix_push_column(cons_jac, 0));

  SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac, 0, 0, 2*data->x[0]));
  SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac, 1, 0, 2*(data->x[0] - 1.)));

  SLEQP_CALL(sleqp_sparse_matrix_push_column(cons_jac, 1));

  SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac, 0, 1, 2*data->x[1]));
  SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac, 1, 1, 2*(data->x[1] - 1.)));

  return SLEQP_OKAY;
}

SLEQP_RETCODE quadconsfunc_hess_prod(SleqpFunc* func,
                                     const double* func_dual,
                                     const SleqpSparseVec* direction,
                                     const SleqpSparseVec* cons_duals,
                                     SleqpSparseVec* result,
                                     void* func_data)
{
  double total_value = 0.;

  if(func_dual)
  {
    total_value += *func_dual;
  }

  for(int k = 0; k < cons_duals->nnz; ++k)
  {
    total_value += cons_duals->data[k];
  }

  total_value *= 2;

  SLEQP_CALL(sleqp_sparse_vector_copy(direction, result));

  SLEQP_CALL(sleqp_sparse_vector_scale(result, total_value));

  return SLEQP_OKAY;
}


void quadconsfunc_setup()
{
  const double inf = sleqp_infinity();

  ASSERT_CALL(sleqp_malloc(&func_data));

  ASSERT_CALL(sleqp_alloc_array(&func_data->x, 2));

  SleqpFuncCallbacks callbacks = {
    .set_value = quadconsfunc_set,
    .func_val = quadconsfunc_val,
    .func_grad = quadconsfunc_grad,
    .cons_val = quadconsfunc_cons_val,
    .cons_jac = quadconsfunc_cons_jac,
    .hess_prod = quadconsfunc_hess_prod,
    .func_free = NULL
  };

  ASSERT_CALL(sleqp_func_create(&quadconsfunc,
                                &callbacks,
                                num_variables,
                                num_constraints,
                                func_data));

  ASSERT_CALL(sleqp_sparse_vector_create(&quadconsfunc_var_lb,
                                         2,
                                         2));

  ASSERT_CALL(sleqp_sparse_vector_push(quadconsfunc_var_lb, 0, 0.));
  ASSERT_CALL(sleqp_sparse_vector_push(quadconsfunc_var_lb, 1, 0.));

  ASSERT_CALL(sleqp_sparse_vector_create(&quadconsfunc_var_ub,
                                         2,
                                         2));

  ASSERT_CALL(sleqp_sparse_vector_push(quadconsfunc_var_ub, 0, 1.));
  ASSERT_CALL(sleqp_sparse_vector_push(quadconsfunc_var_ub, 1, 1.));

  ASSERT_CALL(sleqp_sparse_vector_create(&quadconsfunc_cons_lb,
                                         2,
                                         2));

  ASSERT_CALL(sleqp_sparse_vector_push(quadconsfunc_cons_lb, 0, -inf));
  ASSERT_CALL(sleqp_sparse_vector_push(quadconsfunc_cons_lb, 1, -inf));

  ASSERT_CALL(sleqp_sparse_vector_create(&quadconsfunc_cons_ub,
                                         2,
                                         2));

  ASSERT_CALL(sleqp_sparse_vector_push(quadconsfunc_cons_ub, 0, 1.));
  ASSERT_CALL(sleqp_sparse_vector_push(quadconsfunc_cons_ub, 1, 1.));

  ASSERT_CALL(sleqp_sparse_vector_create(&quadconsfunc_x,
                                         2,
                                         2));

  double val = 0.29289321881345254;

  ASSERT_CALL(sleqp_sparse_vector_push(quadconsfunc_x, 0, val));
  ASSERT_CALL(sleqp_sparse_vector_push(quadconsfunc_x, 1, val));
}

void quadconsfunc_teardown()
{
  ASSERT_CALL(sleqp_sparse_vector_free(&quadconsfunc_x));

  ASSERT_CALL(sleqp_sparse_vector_free(&quadconsfunc_cons_ub));

  ASSERT_CALL(sleqp_sparse_vector_free(&quadconsfunc_cons_lb));

  ASSERT_CALL(sleqp_sparse_vector_free(&quadconsfunc_var_ub));

  ASSERT_CALL(sleqp_sparse_vector_free(&quadconsfunc_var_lb));


  ASSERT_CALL(sleqp_func_release(&quadconsfunc));

  sleqp_free(&func_data->x);

  sleqp_free(&func_data);
}
