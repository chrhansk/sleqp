#include "rosenbrock_lsq_fixture.h"

#include <math.h>

#include "lsq.h"
#include "params.h"


typedef struct RosenbrockData
{
  double a, b;

  double* x;
  double* d;
} RosenbrockData;

RosenbrockData* rosenbrock_func_data;


static inline double sq(double v)
{
  return v*v;
}

static const int num_variables = 2;
static const int num_constraints = 0;
static const int num_residuals = 3;

static SleqpParams* params;

SleqpFunc* rosenbrock_lsq_func;

static SLEQP_RETCODE rosenbrock_lsq_set(SleqpFunc* func,
                                        SleqpSparseVec* x,
                                        SLEQP_VALUE_REASON reason,
                                        bool* reject,
                                        int* func_grad_nnz,
                                        int* cons_val_nnz,
                                        int* cons_jac_nnz,
                                        void* func_data)
{
  *func_grad_nnz = 2;
  *cons_val_nnz = 0;
  *cons_jac_nnz = 0;

  RosenbrockData* data = (RosenbrockData*) func_data;

  SLEQP_CALL(sleqp_sparse_vector_to_raw(x, data->x));

  return SLEQP_OKAY;
}

SLEQP_RETCODE rosenbrock_lsq_eval(SleqpFunc* func,
                                  SleqpSparseVec* residual,
                                  void* func_data)
{
  assert(residual->dim == num_residuals);

  RosenbrockData* data = (RosenbrockData*) func_data;

  const double a = data->a;
  const double b = data->b;

  const double x0 = data->x[0];
  const double x1 = data->x[1];

  SLEQP_CALL(sleqp_sparse_vector_reserve(residual, 2));

  SLEQP_CALL(sleqp_sparse_vector_push(residual, 0, a - x0));

  SLEQP_CALL(sleqp_sparse_vector_push(residual,
                                      1,
                                      sqrt(b)*(x1 - sq(x0))));


  return SLEQP_OKAY;
}

SLEQP_RETCODE rosenbrock_lsq_jac_forward(SleqpFunc* func,
                                         const SleqpSparseVec* forward_direction,
                                         SleqpSparseVec* product,
                                         void* func_data)
{
  assert(forward_direction->dim == num_variables);
  assert(product->dim == num_residuals);

  RosenbrockData* data = (RosenbrockData*) func_data;

  SLEQP_CALL(sleqp_sparse_vector_to_raw(forward_direction, data->d));

  const double b = data->b;

  const double x0 = data->x[0];

  const double d0 = data->d[0];
  const double d1 = data->d[1];

  SLEQP_CALL(sleqp_sparse_vector_reserve(product, 2));

  SLEQP_CALL(sleqp_sparse_vector_push(product,
                                      0,
                                      -1. * d0));

  SLEQP_CALL(sleqp_sparse_vector_push(product,
                                      1,
                                      sqrt(b)* (-2*x0*d0 + d1)));

  return SLEQP_OKAY;
}

SLEQP_RETCODE rosenbrock_lsq_jac_adjoint(SleqpFunc* func,
                                         const SleqpSparseVec* adjoint_direction,
                                         SleqpSparseVec* product,
                                         void* func_data)
{
  assert(adjoint_direction->dim == num_residuals);
  assert(product->dim == num_variables);

  RosenbrockData* data = (RosenbrockData*) func_data;

  SLEQP_CALL(sleqp_sparse_vector_to_raw(adjoint_direction, data->d));

  const double b = data->b;

  const double x0 = data->x[0];

  const double d0 = data->d[0];
  const double d1 = data->d[1];

  SLEQP_CALL(sleqp_sparse_vector_reserve(product, 2));

  SLEQP_CALL(sleqp_sparse_vector_push(product,
                                      0,
                                      -1. * d0 - 2*sqrt(b)*x0*d1));

  SLEQP_CALL(sleqp_sparse_vector_push(product,
                                      1,
                                      sqrt(b)*d1));

  return SLEQP_OKAY;
}

void rosenbrock_lsq_setup()
{
  ASSERT_CALL(sleqp_params_create(&params));

  ASSERT_CALL(sleqp_malloc(&rosenbrock_func_data));

  ASSERT_CALL(sleqp_alloc_array(&rosenbrock_func_data->x, num_variables));
  ASSERT_CALL(sleqp_alloc_array(&rosenbrock_func_data->d, num_residuals));

  rosenbrock_func_data->a = 1.;
  rosenbrock_func_data->b = 100.;

  SleqpLSQCallbacks callbacks = {
    .set_value       = rosenbrock_lsq_set,
    .lsq_residuals   = rosenbrock_lsq_eval,
    .lsq_jac_forward = rosenbrock_lsq_jac_forward,
    .lsq_jac_adjoint = rosenbrock_lsq_jac_adjoint,
    .cons_val        = NULL,
    .cons_jac        = NULL,
    .func_free       = NULL
  };

  ASSERT_CALL(sleqp_lsq_func_create(&rosenbrock_lsq_func,
                                    &callbacks,
                                    num_variables,   // num variables
                                    num_constraints, // num constraints
                                    num_residuals,   // num residuals
                                    0.,              // ML-term
                                    params,
                                    rosenbrock_func_data));
}

void rosenbrock_lsq_teardown()
{
  ASSERT_CALL(sleqp_func_release(&rosenbrock_lsq_func));

  sleqp_free(&rosenbrock_func_data->d);

  sleqp_free(&rosenbrock_func_data->x);

  sleqp_free(&rosenbrock_func_data);

  ASSERT_CALL(sleqp_params_release(&params));
}
