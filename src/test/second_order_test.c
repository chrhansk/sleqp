#include <stdlib.h>
#include <check.h>

#include "sleqp_cmp.h"
#include "sleqp_mem.h"
#include "sleqp_solver.h"

#include "test_common.h"

SleqpSparseVec* var_lb;
SleqpSparseVec* var_ub;
SleqpSparseVec* cons_lb;
SleqpSparseVec* cons_ub;
SleqpSparseVec* x;

const int num_variables = 2;
const int num_constraints = 1;

typedef struct FuncData
{
  double* values;
  double* duals;
  double* direction;

} FuncData;

static double sq(double x)
{
  return x*x;
}

static SLEQP_RETCODE func_set(SleqpSparseVec* x,
                              int num_variables,
                              int* func_grad_nnz,
                              int* cons_val_nnz,
                              int* cons_jac_nnz,
                              void* func_data)
{
  FuncData* data = (FuncData*) func_data;

  SLEQP_CALL(sleqp_sparse_vector_to_raw(x, data->values));

  *func_grad_nnz = num_variables;
  *cons_val_nnz = num_constraints;
  *cons_jac_nnz = num_constraints * num_variables;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE func_eval(int num_variables,
                               SleqpSparseVec* cons_indices,
                               double* func_val,
                               SleqpSparseVec* func_grad,
                               SleqpSparseVec* cons_val,
                               SleqpSparseMatrix* cons_jac,
                               void* func_data)
{
  FuncData* data = (FuncData*) func_data;

  const double x = data->values[0];
  const double y = data->values[1];

  const double xsq = sq(x);
  const double ysq = sq(y);


  if(func_val)
  {
    (*func_val) = 2*(xsq + ysq - 1) - x;
  }

  if(cons_val)
  {
    SLEQP_CALL(sleqp_sparse_vector_push(cons_val,
                                        0,
                                        xsq + ysq - 1));
  }

  if(func_grad)
  {
    SLEQP_CALL(sleqp_sparse_vector_push(func_grad,
                                        0,
                                        4*x - 1));

    SLEQP_CALL(sleqp_sparse_vector_push(func_grad,
                                        1,
                                        4*y));
  }

  if(cons_jac)
  {
    SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac,
                                        0,
                                        0,
                                        2*x));

    SLEQP_CALL(sleqp_sparse_matrix_push_column(cons_jac, 1));

    SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac,
                                        0,
                                        1,
                                        2*y));
  }


  return SLEQP_OKAY;
}

static SLEQP_RETCODE func_hess_prod(int num_variables,
                                    double* func_dual,
                                    SleqpSparseVec* direction,
                                    SleqpSparseVec* cons_duals,
                                    SleqpSparseVec* product,
                                    void* func_data)
{
  FuncData* data = (FuncData*) func_data;

  SLEQP_CALL(sleqp_sparse_vector_to_raw(cons_duals, data->duals));

  SLEQP_CALL(sleqp_sparse_vector_to_raw(direction, data->direction));

  const double d_x = data->direction[0];
  const double d_y = data->direction[1];

  double* duals = data->duals;

  double c_dual = duals[0];

  double f_dual = func_dual ? *func_dual : 0.;

  SLEQP_CALL(sleqp_sparse_vector_reserve(product, num_variables));

  {
    double v = (4*f_dual + 2*c_dual)*d_x;

    SLEQP_CALL(sleqp_sparse_vector_push(product, 0, v));
  }

  {
    double v = (4*f_dual + 2*c_dual)*d_y;

    SLEQP_CALL(sleqp_sparse_vector_push(product, 1, v));
  }

  return SLEQP_OKAY;
}

FuncData* func_data;
SleqpFunc* func;

SleqpParams* params;
SleqpProblem* problem;

void second_order_setup()
{
  const double inf = sleqp_infinity();

  ASSERT_CALL(sleqp_sparse_vector_create(&var_lb,
                                         num_variables,
                                         num_variables));

  ASSERT_CALL(sleqp_sparse_vector_create(&var_ub,
                                         num_variables,
                                         num_variables));

  for(int i = 0; i < num_variables; ++i)
  {
    ASSERT_CALL(sleqp_sparse_vector_push(var_lb, i, -inf));
    ASSERT_CALL(sleqp_sparse_vector_push(var_ub, i, inf));
  }

  ASSERT_CALL(sleqp_sparse_vector_create(&cons_lb,
                                         num_constraints,
                                         num_constraints));

  ASSERT_CALL(sleqp_sparse_vector_create(&cons_ub,
                                         num_constraints,
                                         num_constraints));

  ASSERT_CALL(sleqp_sparse_vector_push(cons_lb, 0, 0.));
  ASSERT_CALL(sleqp_sparse_vector_push(cons_ub, 0, 0.));

  ASSERT_CALL(sleqp_sparse_vector_create(&x,
                                         num_variables,
                                         num_variables));

  ASSERT_CALL(sleqp_sparse_vector_push(x, 0, 0.707106781186548));
  ASSERT_CALL(sleqp_sparse_vector_push(x, 1, 0.707106781186548));

  ASSERT_CALL(sleqp_malloc(&func_data));

  ASSERT_CALL(sleqp_calloc(&func_data->values, num_variables));
  ASSERT_CALL(sleqp_calloc(&func_data->duals, num_constraints));
  ASSERT_CALL(sleqp_calloc(&func_data->direction, num_variables));

  SleqpFuncCallbacks callbacks = {
    .set_value = func_set,
    .func_eval = func_eval,
    .hess_prod = func_hess_prod,
    .func_free = NULL
  };

  ASSERT_CALL(sleqp_func_create(&func,
                                &callbacks,
                                num_variables,
                                func_data));

  ASSERT_CALL(sleqp_params_create(&params));

  ASSERT_CALL(sleqp_problem_create(&problem,
                                   func,
                                   params,
                                   var_lb,
                                   var_ub,
                                   cons_lb,
                                   cons_ub));
}

void second_order_teardown()
{
  ASSERT_CALL(sleqp_problem_free(&problem));

  ASSERT_CALL(sleqp_params_free(&params));

  ASSERT_CALL(sleqp_func_free(&func));

  sleqp_free(&func_data->direction);
  sleqp_free(&func_data->duals);
  sleqp_free(&func_data->values);
  sleqp_free(&func_data);

  ASSERT_CALL(sleqp_sparse_vector_free(&x));
  ASSERT_CALL(sleqp_sparse_vector_free(&cons_ub));
  ASSERT_CALL(sleqp_sparse_vector_free(&cons_lb));
  ASSERT_CALL(sleqp_sparse_vector_free(&var_ub));
  ASSERT_CALL(sleqp_sparse_vector_free(&var_lb));
}

START_TEST(test_second_order_solve)
{
  SleqpSparseVec* expected_solution;

  SleqpParams* params;
  SleqpOptions* options;
  SleqpSolver* solver;

  ASSERT_CALL(sleqp_sparse_vector_create(&expected_solution, 2, 2));

  ASSERT_CALL(sleqp_sparse_vector_push(expected_solution, 0, 1.));

  ASSERT_CALL(sleqp_params_create(&params));

  ASSERT_CALL(sleqp_options_create(&options));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  params,
                                  options,
                                  x,
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

  ASSERT_CALL(sleqp_options_free(&options));

  ASSERT_CALL(sleqp_params_free(&params));

  ASSERT_CALL(sleqp_sparse_vector_free(&expected_solution));
}
END_TEST

Suite* second_order_test_suite()
{
  Suite *suite;
  TCase *tc_cons;

  suite = suite_create("Constrained tests");

  tc_cons = tcase_create("Constrained solution test");

  tcase_add_checked_fixture(tc_cons,
                            second_order_setup,
                            second_order_teardown);

  tcase_add_test(tc_cons, test_second_order_solve);
  suite_add_tcase(suite, tc_cons);

  return suite;
}


int main()
{
  int num_fails;
  Suite* suite;
  SRunner* srunner;

  suite = second_order_test_suite();
  srunner = srunner_create(suite);

  srunner_set_fork_status(srunner, CK_NOFORK);
  srunner_run_all(srunner, CK_NORMAL);

  num_fails = srunner_ntests_failed(srunner);

  srunner_free(srunner);

  return (num_fails > 0) ? EXIT_FAILURE : EXIT_SUCCESS;
}
