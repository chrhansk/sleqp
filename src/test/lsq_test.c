#include <stdlib.h>
#include <check.h>

#include "sleqp.h"
#include "sleqp_cauchy.h"
#include "sleqp_cmp.h"
#include "sleqp_lsq.h"
#include "sleqp_mem.h"

#include "lp/sleqp_lpi.h"

#include "test_common.h"

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
static const int num_residuals = 3;

SleqpFunc* rosenbrock_lsq_func;

SleqpSparseVec* rosenbrock_var_lb;
SleqpSparseVec* rosenbrock_var_ub;
SleqpSparseVec* rosenbrock_cons_lb;
SleqpSparseVec* rosenbrock_cons_ub;
SleqpSparseVec* rosenbrock_x;

static SLEQP_RETCODE rosenbrock_lsq_set(SleqpSparseVec* x,
                                        int num_variables,
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

SLEQP_RETCODE rosenbrock_lsq_eval(int num_variables,
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

SLEQP_RETCODE rosenbrock_lsq_jac_forward(int num_variables,
                                         SleqpSparseVec* forward_direction,
                                         SleqpSparseVec* product,
                                         void* func_data)
{
  assert(forward_direction->dim == num_variables);
  assert(product->dim == num_residuals);

  RosenbrockData* data = (RosenbrockData*) func_data;

  SLEQP_CALL(sleqp_sparse_vector_to_raw(forward_direction, data->d));

  const double a = data->a;
  const double b = data->b;

  const double x0 = data->x[0];
  const double x1 = data->x[1];

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

SLEQP_RETCODE rosenbrock_lsq_jac_adjoint(int num_variables,
                                         SleqpSparseVec* adjoint_direction,
                                         SleqpSparseVec* product,
                                         void* func_data)
{
  assert(adjoint_direction->dim == num_residuals);
  assert(product->dim == num_variables);

  RosenbrockData* data = (RosenbrockData*) func_data;

  SLEQP_CALL(sleqp_sparse_vector_to_raw(adjoint_direction, data->d));

  const double a = data->a;
  const double b = data->b;

  const double x0 = data->x[0];
  const double x1 = data->x[1];

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

void rosenbrock_setup()
{
  const double inf = sleqp_infinity();

  ASSERT_CALL(sleqp_malloc(&rosenbrock_func_data));

  ASSERT_CALL(sleqp_calloc(&rosenbrock_func_data->x, num_variables));
  ASSERT_CALL(sleqp_calloc(&rosenbrock_func_data->d, num_variables));

  rosenbrock_func_data->a = 1.;
  rosenbrock_func_data->b = 100.;

  SleqpLSQCallbacks callbacks = {
    .set_value = rosenbrock_lsq_set,
    .lsq_eval = rosenbrock_lsq_eval,
    .lsq_jac_forward = rosenbrock_lsq_jac_forward,
    .lsq_jac_adjoint = rosenbrock_lsq_jac_adjoint,
    .eval_additional = NULL,
    .hess_prod_additional = NULL,
    .func_free = NULL
  };

  ASSERT_CALL(sleqp_lsq_func_create(&rosenbrock_lsq_func,
                                    &callbacks,
                                    num_variables, // num variables
                                    num_residuals, // num residuals
                                    0.,            // ML-term
                                    rosenbrock_func_data));

  ASSERT_CALL(sleqp_sparse_vector_create(&rosenbrock_var_lb,
                                         num_variables,
                                         num_variables));

  ASSERT_CALL(sleqp_sparse_vector_push(rosenbrock_var_lb, 0, -inf));
  ASSERT_CALL(sleqp_sparse_vector_push(rosenbrock_var_lb, 1, -inf));

  ASSERT_CALL(sleqp_sparse_vector_create(&rosenbrock_var_ub,
                                         num_variables,
                                         num_variables));

  ASSERT_CALL(sleqp_sparse_vector_push(rosenbrock_var_ub, 0, inf));
  ASSERT_CALL(sleqp_sparse_vector_push(rosenbrock_var_ub, 1, inf));

  ASSERT_CALL(sleqp_sparse_vector_create(&rosenbrock_cons_lb,
                                         0,
                                         0));

  ASSERT_CALL(sleqp_sparse_vector_create(&rosenbrock_cons_ub,
                                         0,
                                         0));

  ASSERT_CALL(sleqp_sparse_vector_create(&rosenbrock_x,
                                         num_variables,
                                         num_variables));

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



  ASSERT_CALL(sleqp_func_free(&rosenbrock_lsq_func));

  sleqp_free(&rosenbrock_func_data->d);

  sleqp_free(&rosenbrock_func_data->x);

  sleqp_free(&rosenbrock_func_data);
}

START_TEST(test_unconstrained_solve)
{
  SleqpSparseVec* expected_solution;

  SleqpParams* params;
  SleqpOptions* options;
  SleqpProblem* problem;
  SleqpSolver* solver;


  ASSERT_CALL(sleqp_sparse_vector_create(&expected_solution,
                                         num_variables,
                                         num_variables));

  ASSERT_CALL(sleqp_sparse_vector_push(expected_solution, 0, 1.));
  ASSERT_CALL(sleqp_sparse_vector_push(expected_solution, 1, 1.));

  ASSERT_CALL(sleqp_params_create(&params));

  ASSERT_CALL(sleqp_params_set_optimality_tolerance(params, 1e-8));

  ASSERT_CALL(sleqp_options_create(&options));

  ASSERT_CALL(sleqp_options_set_deriv_check(options, SLEQP_DERIV_CHECK_FIRST));

  ASSERT_CALL(sleqp_problem_create(&problem,
                                   rosenbrock_lsq_func,
                                   params,
                                   rosenbrock_var_lb,
                                   rosenbrock_var_ub,
                                   rosenbrock_cons_lb,
                                   rosenbrock_cons_ub));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  params,
                                  options,
                                  rosenbrock_x,
                                  NULL));

  // 100 iterations should be plenty...
  ASSERT_CALL(sleqp_solver_solve(solver, 100, -1));

  SleqpIterate* solution_iterate;

  ASSERT_CALL(sleqp_solver_get_solution(solver,
                                        &solution_iterate));

  ck_assert_int_eq(sleqp_solver_get_status(solver), SLEQP_OPTIMAL);

  SleqpSparseVec* actual_solution = solution_iterate->primal;

  ck_assert(sleqp_sparse_vector_eq(actual_solution,
                                   expected_solution,
                                   1e-6));

  ASSERT_CALL(sleqp_solver_free(&solver));

  ASSERT_CALL(sleqp_problem_free(&problem));

  ASSERT_CALL(sleqp_options_free(&options));

  ASSERT_CALL(sleqp_params_free(&params));

  ASSERT_CALL(sleqp_sparse_vector_free(&expected_solution));
}
END_TEST

Suite* lsq_test_suite()
{
  Suite *suite;
  TCase *tc_uncons;

  suite = suite_create("LSQ tests");

  tc_uncons = tcase_create("LSQ solution test");

  tcase_add_checked_fixture(tc_uncons,
                            rosenbrock_setup,
                            rosenbrock_teardown);

  tcase_add_test(tc_uncons, test_unconstrained_solve);
  suite_add_tcase(suite, tc_uncons);

  return suite;
}


int main()
{
  int num_fails;
  Suite* suite;
  SRunner* srunner;

  suite = lsq_test_suite();
  srunner = srunner_create(suite);

  srunner_set_fork_status(srunner, CK_NOFORK);
  srunner_run_all(srunner, CK_NORMAL);

  num_fails = srunner_ntests_failed(srunner);

  srunner_free(srunner);

  return (num_fails > 0) ? EXIT_FAILURE : EXIT_SUCCESS;
}
