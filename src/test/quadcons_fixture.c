#include "quadcons_fixture.h"

SleqpFunc* quadconsfunc;

SleqpSparseVec* quadconsfunc_var_lb;
SleqpSparseVec* quadconsfunc_var_ub;
SleqpSparseVec* quadconsfunc_cons_lb;
SleqpSparseVec* quadconsfunc_cons_ub;
SleqpSparseVec* quadconsfunc_x;

typedef struct SquareFuncData
{
  double* x;
} SquareFuncData;

SquareFuncData* func_data;

static inline double square(double v)
{
  return v*v;
}

SLEQP_RETCODE quadconsfunc_set(SleqpSparseVec* x,
                               int num_variables,
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

SLEQP_RETCODE quadconsfunc_eval(int num_variables,
                                int* indices,
                                double* func_val,
                                SleqpSparseVec* func_grad,
                                SleqpSparseVec* cons_val,
                                SleqpSparseMatrix* cons_jac,
                                void* func_data)
{
  SquareFuncData* data = (SquareFuncData*) func_data;

  if(func_val)
  {
    *func_val = square(data->x[0]) + square(data->x[1]);
  }

  if(func_grad)
  {
    assert(func_grad->dim == 2);
    assert(func_grad->nnz_max >= 2);

    func_grad->nnz = 0;

    SLEQP_CALL(sleqp_sparse_vector_push(func_grad,
                                        0,
                                        2.*data->x[0]));

    SLEQP_CALL(sleqp_sparse_vector_push(func_grad,
                                        1,
                                        2.*data->x[1]));

  }

  if(cons_val)
  {
    assert(cons_val->dim == 2);
    assert(cons_val->nnz_max >= 2);

    cons_val->nnz = 0;

    SLEQP_CALL(sleqp_sparse_vector_push(cons_val,
                                        0,
                                        square(data->x[0]) + square(data->x[1])));

    SLEQP_CALL(sleqp_sparse_vector_push(cons_val,
                                        1,
                                        (square(1 - data->x[0])) + square(1 - data->x[1])));

  }

  if(cons_jac)
  {
    assert(cons_jac->num_rows == 2);
    assert(cons_jac->num_cols == 2);
    assert(cons_jac->nnz_max >= 4);

    SLEQP_CALL(sleqp_sparse_matrix_add_column(cons_jac, 0));

    SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac, 0, 0, 2*data->x[0]));
    SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac, 1, 0, 2*(data->x[0] - 1.)));

    SLEQP_CALL(sleqp_sparse_matrix_add_column(cons_jac, 1));

    SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac, 0, 1, 2*data->x[1]));
    SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac, 1, 1, 2*(data->x[1] - 1.)));

    SLEQP_CALL(sleqp_sparse_matrix_fprintf(cons_jac,
                                           stdout));

  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE quadconsfunc_eval_bilinear(int num_variables,
                                         double* func_dual,
                                         SleqpSparseVec* direction,
                                         SleqpSparseVec* cons_duals,
                                         double* bilinear_prod,
                                         void* func_data)
{
  return SLEQP_OKAY;
}


void quadconsfunc_setup()
{
  const double inf = sleqp_infinity();

  ASSERT_CALL(sleqp_malloc(&func_data));

  ASSERT_CALL(sleqp_calloc(&func_data->x, 2));

  ASSERT_CALL(sleqp_func_create(&quadconsfunc,
                                quadconsfunc_set,
                                quadconsfunc_eval,
                                quadconsfunc_eval_bilinear,
                                2,
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


  ASSERT_CALL(sleqp_func_free(&quadconsfunc));

  sleqp_free(&func_data->x);

  sleqp_free(&func_data);
}
