#include "rosenbrock_fixture.h"

const int num_variables = 2;
const int num_constraints = 0;

typedef struct RosenbrockData
{
  double a, b;

  double* x;
} RosenbrockData;


static inline double sq(double v)
{
  return v*v;
}


RosenbrockData* func_data;

SleqpFunc* rosenbrock_func;

SleqpSparseVec* rosenbrock_var_lb;
SleqpSparseVec* rosenbrock_var_ub;
SleqpSparseVec* rosenbrock_cons_lb;
SleqpSparseVec* rosenbrock_cons_ub;
SleqpSparseVec* rosenbrock_x;


static SLEQP_RETCODE rosenbrock_set(SleqpFunc* func,
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

  RosenbrockData* data = (RosenbrockData*) func_data;

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

static SLEQP_RETCODE rosenbrock_eval(SleqpFunc* func,
                                     const SleqpSparseVec* cons_indices,
                                     double* func_val,
                                     SleqpSparseVec* func_grad,
                                     SleqpSparseVec* cons_val,
                                     SleqpSparseMatrix* cons_jac,
                                     void* func_data)
{
  RosenbrockData* data = (RosenbrockData*) func_data;

  double x = data->x[0];
  double y = data->x[1];

  double a = data->a;
  double b = data->b;

  double xsq = sq(x);

  if(func_val)
  {
    *func_val = sq(a - x) + b * sq(y - xsq);
  }

  if(func_grad)
  {
    assert(func_grad->nnz == 0);
    assert(func_grad->dim == 2);

    double gradx = (4.*b*x*(xsq - y)) + 2.*x - 2.*a;

    double grady = -2.*b*(xsq - y);

    SLEQP_CALL(sleqp_sparse_vector_push(func_grad,
                                        0,
                                        gradx));

    SLEQP_CALL(sleqp_sparse_vector_push(func_grad,
                                        1,
                                        grady));

  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE rosenbrock_hess_prod(SleqpFunc* func,
                                          const double* func_dual,
                                          const SleqpSparseVec* direction,
                                          const SleqpSparseVec* cons_duals,
                                          SleqpSparseVec* product,
                                          void* func_data)
{
  RosenbrockData* data = (RosenbrockData*) func_data;

  double x = data->x[0];
  double y = data->x[1];

  double b = data->b;

  double xsq = sq(x);

  double d[2];

  SLEQP_CALL(sleqp_sparse_vector_to_raw(direction, d));

  if(func_dual)
  {
    SLEQP_CALL(sleqp_sparse_vector_reserve(product, 2));

    SLEQP_CALL(sleqp_sparse_vector_push(product,
                                        0,
                                        (8.*b*xsq + 4.*b*(xsq - y) + 2.)*d[0]
                                        - (4.*b*x)*d[1]));

    SLEQP_CALL(sleqp_sparse_vector_push(product,
                                        1,
                                        (-4.*b*x)*d[0] + (2.*b)*d[1]));


  }

  return SLEQP_OKAY;
}


void rosenbrock_setup()
{
  const double inf = sleqp_infinity();

  ASSERT_CALL(sleqp_malloc(&func_data));

  func_data->a = 1.;
  func_data->b = 100.;

  ASSERT_CALL(sleqp_calloc(&func_data->x, 2));

    SleqpFuncCallbacks callbacks = {
    .set_value = rosenbrock_set,
    .func_eval = rosenbrock_eval,
    .hess_prod = rosenbrock_hess_prod,
    .func_free = NULL
  };

  ASSERT_CALL(sleqp_func_create(&rosenbrock_func,
                                &callbacks,
                                num_variables,
                                num_constraints,
                                func_data));

  ASSERT_CALL(sleqp_sparse_vector_create(&rosenbrock_var_lb,
                                         2,
                                         2));

  ASSERT_CALL(sleqp_sparse_vector_push(rosenbrock_var_lb, 0, -inf));
  ASSERT_CALL(sleqp_sparse_vector_push(rosenbrock_var_lb, 1, -inf));

  ASSERT_CALL(sleqp_sparse_vector_create(&rosenbrock_var_ub,
                                         2,
                                         2));

  ASSERT_CALL(sleqp_sparse_vector_push(rosenbrock_var_ub, 0, inf));
  ASSERT_CALL(sleqp_sparse_vector_push(rosenbrock_var_ub, 1, inf));

  ASSERT_CALL(sleqp_sparse_vector_create(&rosenbrock_cons_lb,
                                         0,
                                         0));

  ASSERT_CALL(sleqp_sparse_vector_create(&rosenbrock_cons_ub,
                                         0,
                                         0));

  ASSERT_CALL(sleqp_sparse_vector_create(&rosenbrock_x,
                                         2,
                                         2));

  ASSERT_CALL(sleqp_sparse_vector_push(rosenbrock_x, 0, 0.));
  ASSERT_CALL(sleqp_sparse_vector_push(rosenbrock_x, 1, 0.));
}

void rosenbrock_teardown()
{
  ASSERT_CALL(sleqp_sparse_vector_free(&rosenbrock_x));

  ASSERT_CALL(sleqp_sparse_vector_free(&rosenbrock_cons_ub));

  ASSERT_CALL(sleqp_sparse_vector_free(&rosenbrock_cons_lb));

  ASSERT_CALL(sleqp_sparse_vector_free(&rosenbrock_var_ub));

  ASSERT_CALL(sleqp_sparse_vector_free(&rosenbrock_var_lb));



  ASSERT_CALL(sleqp_func_release(&rosenbrock_func));

  sleqp_free(&func_data->x);

  sleqp_free(&func_data);
}
