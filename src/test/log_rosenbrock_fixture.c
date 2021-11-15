#include "log_rosenbrock_fixture.h"

#include <math.h>

const int log_rosenbrock_num_variables   = 2;
const int log_rosenbrock_num_constraints = 0;

SleqpFunc* log_rosenbrock_func;

SleqpSparseVec* log_rosenbrock_var_lb;
SleqpSparseVec* log_rosenbrock_var_ub;
SleqpSparseVec* log_rosenbrock_cons_lb;
SleqpSparseVec* log_rosenbrock_cons_ub;
SleqpSparseVec* log_rosenbrock_initial;
SleqpSparseVec* log_rosenbrock_optimal;

static double x;
static double y;

static double inner;

static double
sq(double x)
{
  return x * x;
}

static double
eval()
{
  return 1. + 10000. * sq(y - sq(x)) + sq(1 - x);
}

static SLEQP_RETCODE
log_rosenbrock_set(SleqpFunc* func,
                   SleqpSparseVec* v,
                   SLEQP_VALUE_REASON reason,
                   bool* reject,
                   int* func_grad_nnz,
                   int* cons_val_nnz,
                   int* cons_jac_nnz,
                   void* func_data)
{
  x     = sleqp_sparse_vector_value_at(v, 0);
  y     = sleqp_sparse_vector_value_at(v, 1);
  inner = eval();

  *func_grad_nnz = 2;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
log_rosenbrock_val(SleqpFunc* func, double* func_val, void* func_data)
{
  *func_val = log(inner);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
log_rosenbrock_grad(SleqpFunc* func, SleqpSparseVec* func_grad, void* func_data)
{
  const double dx = 2. * (20000. * (x * sq(x) - x * y) + x - 1) / inner;

  SLEQP_CALL(sleqp_sparse_vector_push(func_grad, 0, dx));

  const double dy = 20000. * (y - sq(x)) / inner;

  SLEQP_CALL(sleqp_sparse_vector_push(func_grad, 1, dy));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
log_rosenbrock_hess_prod(SleqpFunc* func,
                         const double* func_dual,
                         const SleqpSparseVec* direction,
                         const SleqpSparseVec* cons_duals,
                         SleqpSparseVec* product,
                         void* func_data)
{
  SLEQP_CALL(
    sleqp_sparse_vector_reserve(product, log_rosenbrock_num_variables));

  const double dx = sleqp_sparse_vector_value_at(direction, 0);
  const double dy = sleqp_sparse_vector_value_at(direction, 1);

  const double hxx
    = (-40000. * (y - sq(x)) + 80000. * sq(x) + 2) / inner
      - sq(-40000. * x * (y - sq(x)) - 2 * sq(1 - x)) / sq(inner);

  const double hxy
    = (-40000. * x) / inner
      - (20000. * (y - sq(x)) * (-40000. * x * (y - sq(x)) - 2 * (1 - x)))
          / sq(inner);

  const double hyy = 20000. / inner + sq(20000.) * sq(y - sq(x)) / sq(inner);

  if (func_dual)
  {
    SLEQP_CALL(
      sleqp_sparse_vector_push(product, 0, (*func_dual) * hxx * dx + hxy * dy));

    SLEQP_CALL(
      sleqp_sparse_vector_push(product, 1, (*func_dual) * hxy * dx + hyy * dy));
  }

  return SLEQP_OKAY;
}

void
log_rosenbrock_setup()
{
  const double inf = sleqp_infinity();

  SleqpFuncCallbacks callbacks = {.set_value = log_rosenbrock_set,
                                  .func_val  = log_rosenbrock_val,
                                  .func_grad = log_rosenbrock_grad,
                                  .cons_val  = NULL,
                                  .cons_jac  = NULL,
                                  .hess_prod = log_rosenbrock_hess_prod,
                                  .func_free = NULL};

  ASSERT_CALL(sleqp_func_create(&log_rosenbrock_func,
                                &callbacks,
                                log_rosenbrock_num_variables,
                                log_rosenbrock_num_constraints,
                                NULL));

  ASSERT_CALL(sleqp_sparse_vector_create(&log_rosenbrock_var_lb, 2, 2));

  ASSERT_CALL(sleqp_sparse_vector_push(log_rosenbrock_var_lb, 0, -inf));
  ASSERT_CALL(sleqp_sparse_vector_push(log_rosenbrock_var_lb, 1, -inf));

  ASSERT_CALL(sleqp_sparse_vector_create(&log_rosenbrock_var_ub, 2, 2));

  ASSERT_CALL(sleqp_sparse_vector_push(log_rosenbrock_var_ub, 0, inf));
  ASSERT_CALL(sleqp_sparse_vector_push(log_rosenbrock_var_ub, 1, inf));

  ASSERT_CALL(sleqp_sparse_vector_create_empty(&log_rosenbrock_cons_lb, 0));

  ASSERT_CALL(sleqp_sparse_vector_create_empty(&log_rosenbrock_cons_ub, 0));

  ASSERT_CALL(sleqp_sparse_vector_create_empty(&log_rosenbrock_initial, 2));

  ASSERT_CALL(sleqp_sparse_vector_create_full(&log_rosenbrock_optimal, 2));

  ASSERT_CALL(sleqp_sparse_vector_push(log_rosenbrock_optimal, 0, 1.));
  ASSERT_CALL(sleqp_sparse_vector_push(log_rosenbrock_optimal, 1, 1.));
}

void
log_rosenbrock_teardown()
{
  ASSERT_CALL(sleqp_sparse_vector_free(&log_rosenbrock_optimal));

  ASSERT_CALL(sleqp_sparse_vector_free(&log_rosenbrock_initial));

  ASSERT_CALL(sleqp_sparse_vector_free(&log_rosenbrock_cons_ub));

  ASSERT_CALL(sleqp_sparse_vector_free(&log_rosenbrock_cons_lb));

  ASSERT_CALL(sleqp_sparse_vector_free(&log_rosenbrock_var_ub));

  ASSERT_CALL(sleqp_sparse_vector_free(&log_rosenbrock_var_lb));

  ASSERT_CALL(sleqp_func_release(&log_rosenbrock_func));
}
