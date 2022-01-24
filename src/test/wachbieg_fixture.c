#include "wachbieg_fixture.h"

#include "test_common.h"

#include "cmp.h"
#include "mem.h"

#include "sparse/sparse_matrix.h"

const int wachbieg_num_variables   = 3;
const int wachbieg_num_constraints = 2;

const double a = -1.;
const double b = .5;

SleqpFunc* wachbieg_func;

SleqpSparseVec* wachbieg_var_lb;
SleqpSparseVec* wachbieg_var_ub;
SleqpSparseVec* wachbieg_cons_lb;
SleqpSparseVec* wachbieg_cons_ub;
SleqpSparseVec* wachbieg_initial;
SleqpSparseVec* wachbieg_optimal;

SleqpFunc* wachbieg_func;

typedef struct
{
  double x[3];
} WachBiegData;

static inline double
sq(double v)
{
  return v * v;
}

WachBiegData wachbieg_data;

const int wachbieg_jac_nnz = 4;

static SLEQP_RETCODE
wachbieg_set(SleqpFunc* func,
             SleqpSparseVec* x,
             SLEQP_VALUE_REASON reason,
             bool* reject,
             int* obj_grad_nnz,
             int* cons_val_nnz,
             int* cons_jac_nnz,
             void* func_data)
{
  *obj_grad_nnz = wachbieg_num_variables;
  *cons_val_nnz = wachbieg_num_constraints;
  *cons_jac_nnz = wachbieg_jac_nnz;

  wachbieg_data.x[0] = 0;
  wachbieg_data.x[1] = 0;
  wachbieg_data.x[2] = 0;

  int k_x = 0;

  while (k_x < x->nnz)
  {
    wachbieg_data.x[x->indices[k_x]] = x->data[k_x];

    ++k_x;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
wachbieg_obj_val(SleqpFunc* func, double* obj_val, void* func_data)
{
  const double x0 = wachbieg_data.x[0];

  *obj_val = x0;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
wachbieg_obj_grad(SleqpFunc* func, SleqpSparseVec* obj_grad, void* func_data)
{
  SLEQP_CALL(sleqp_sparse_vector_push(obj_grad, 0, 1.));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
wachbieg_cons_val(SleqpFunc* func, SleqpSparseVec* cons_val, void* func_data)
{
  const double x0 = wachbieg_data.x[0];
  const double x1 = wachbieg_data.x[1];
  const double x2 = wachbieg_data.x[2];

  SLEQP_CALL(sleqp_sparse_vector_push(cons_val, 0, sq(x0) - x1 + a));

  SLEQP_CALL(sleqp_sparse_vector_push(cons_val, 1, x0 - x2 - b));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
wachbieg_cons_jac(SleqpFunc* func, SleqpSparseMatrix* cons_jac, void* func_data)
{
  const double x0 = wachbieg_data.x[0];

  SLEQP_CALL(sleqp_sparse_matrix_clear(cons_jac));

  SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac, 0, 0, 2 * x0));

  SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac, 1, 0, 1.));

  SLEQP_CALL(sleqp_sparse_matrix_push_column(cons_jac, 1));

  SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac, 0, 1, -1.));

  SLEQP_CALL(sleqp_sparse_matrix_push_column(cons_jac, 2));

  SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac, 1, 2, -1.));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
wachbieg_hess_prod(SleqpFunc* func,
                   const double* obj_dual,
                   const SleqpSparseVec* direction,
                   const SleqpSparseVec* cons_duals,
                   SleqpSparseVec* product,
                   void* func_data)
{
  const double d0 = sleqp_sparse_vector_value_at(direction, 0);
  const double m  = sleqp_sparse_vector_value_at(cons_duals, 0);

  SLEQP_CALL(sleqp_sparse_vector_reserve(product, 1));

  SLEQP_CALL(sleqp_sparse_vector_push(product, 0, 2. * d0 * m));

  return SLEQP_OKAY;
}

void
wachbieg_setup()
{
  const double inf = sleqp_infinity();

  SleqpFuncCallbacks callbacks = {.set_value = wachbieg_set,
                                  .obj_val   = wachbieg_obj_val,
                                  .obj_grad  = wachbieg_obj_grad,
                                  .cons_val  = wachbieg_cons_val,
                                  .cons_jac  = wachbieg_cons_jac,
                                  .hess_prod = wachbieg_hess_prod,
                                  .func_free = NULL};

  ASSERT_CALL(sleqp_func_create(&wachbieg_func,
                                &callbacks,
                                wachbieg_num_variables,
                                wachbieg_num_constraints,
                                NULL));

  ASSERT_CALL(
    sleqp_sparse_vector_create_full(&wachbieg_var_lb, wachbieg_num_variables));

  {
    double values[3] = {-inf, 0., 0.};

    ASSERT_CALL(sleqp_sparse_vector_from_raw(wachbieg_var_lb,
                                             values,
                                             wachbieg_num_variables,
                                             0.));
  }

  ASSERT_CALL(
    sleqp_sparse_vector_create_full(&wachbieg_var_ub, wachbieg_num_variables));

  {
    double values[3] = {inf, inf, inf};

    ASSERT_CALL(sleqp_sparse_vector_from_raw(wachbieg_var_ub,
                                             values,
                                             wachbieg_num_variables,
                                             0.));
  }

  ASSERT_CALL(sleqp_sparse_vector_create_empty(&wachbieg_cons_lb,
                                               wachbieg_num_constraints));

  ASSERT_CALL(sleqp_sparse_vector_create_empty(&wachbieg_cons_ub,
                                               wachbieg_num_constraints));

  ASSERT_CALL(
    sleqp_sparse_vector_create_full(&wachbieg_initial, wachbieg_num_variables));

  {
    double values[3] = {-2., 1., 1.};

    ASSERT_CALL(sleqp_sparse_vector_from_raw(wachbieg_initial,
                                             values,
                                             wachbieg_num_variables,
                                             0.));
  }

  ASSERT_CALL(
    sleqp_sparse_vector_create_full(&wachbieg_optimal, wachbieg_num_variables));

  // Note: This is not unique...
  // General solutions: [x1, x2, x3] with
  // x1 = x3 + b, x2 = (x3 + b)^2 + a, x3 >= 0 arbitrary
  {
    double values[3] = {1., 0., .5};

    ASSERT_CALL(sleqp_sparse_vector_from_raw(wachbieg_optimal,
                                             values,
                                             wachbieg_num_variables,
                                             0.));
  }
}

void
wachbieg_teardown()
{
  ASSERT_CALL(sleqp_sparse_vector_free(&wachbieg_optimal));
  ASSERT_CALL(sleqp_sparse_vector_free(&wachbieg_initial));

  ASSERT_CALL(sleqp_sparse_vector_free(&wachbieg_cons_ub));
  ASSERT_CALL(sleqp_sparse_vector_free(&wachbieg_cons_lb));

  ASSERT_CALL(sleqp_sparse_vector_free(&wachbieg_var_ub));
  ASSERT_CALL(sleqp_sparse_vector_free(&wachbieg_var_lb));

  ASSERT_CALL(sleqp_func_release(&wachbieg_func));
}
