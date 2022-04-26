#include <check.h>
#include <stdlib.h>

#include "cmp.h"
#include "mem.h"
#include "solver.h"

#include "test_common.h"

SleqpVec* var_lb;
SleqpVec* var_ub;
SleqpVec* cons_lb;
SleqpVec* cons_ub;
SleqpVec* x;

static const int num_variables   = 2;
static const int num_constraints = 1;

typedef struct FuncData
{
  double* values;
  double* duals;
  double* direction;

} FuncData;

static double
sq(double x)
{
  return x * x;
}

static SLEQP_RETCODE
func_set(SleqpFunc* func,
         SleqpVec* x,
         SLEQP_VALUE_REASON reason,
         bool* reject,
         void* func_data)
{
  FuncData* data = (FuncData*)func_data;

  SLEQP_CALL(sleqp_vec_to_raw(x, data->values));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
func_obj_val(SleqpFunc* func, double* obj_val, void* func_data)
{
  FuncData* data = (FuncData*)func_data;

  const double x = data->values[0];
  const double y = data->values[1];

  const double xsq = sq(x);
  const double ysq = sq(y);

  (*obj_val) = 2 * (xsq + ysq - 1) - x;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
func_obj_grad(SleqpFunc* func, SleqpVec* obj_grad, void* func_data)
{
  FuncData* data = (FuncData*)func_data;

  const double x = data->values[0];
  const double y = data->values[1];

  SLEQP_CALL(sleqp_vec_push(obj_grad, 0, 4 * x - 1));

  SLEQP_CALL(sleqp_vec_push(obj_grad, 1, 4 * y));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
func_cons_val(SleqpFunc* func, SleqpVec* cons_val, void* func_data)
{
  FuncData* data = (FuncData*)func_data;

  const double x = data->values[0];
  const double y = data->values[1];

  const double xsq = sq(x);
  const double ysq = sq(y);

  SLEQP_CALL(sleqp_vec_push(cons_val, 0, xsq + ysq - 1));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
func_cons_jac(SleqpFunc* func, SleqpSparseMatrix* cons_jac, void* func_data)
{
  FuncData* data = (FuncData*)func_data;

  const double x = data->values[0];
  const double y = data->values[1];

  SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac, 0, 0, 2 * x));

  SLEQP_CALL(sleqp_sparse_matrix_push_column(cons_jac, 1));

  SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac, 0, 1, 2 * y));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
func_hess_prod(SleqpFunc* func,
               const double* obj_dual,
               const SleqpVec* direction,
               const SleqpVec* cons_duals,
               SleqpVec* product,
               void* func_data)
{
  FuncData* data = (FuncData*)func_data;

  SLEQP_CALL(sleqp_vec_to_raw(cons_duals, data->duals));

  SLEQP_CALL(sleqp_vec_to_raw(direction, data->direction));

  const double d_x = data->direction[0];
  const double d_y = data->direction[1];

  double* duals = data->duals;

  double c_dual = duals[0];

  double o_dual = obj_dual ? *obj_dual : 0.;

  SLEQP_CALL(sleqp_vec_reserve(product, num_variables));

  {
    double v = (4 * o_dual + 2 * c_dual) * d_x;

    SLEQP_CALL(sleqp_vec_push(product, 0, v));
  }

  {
    double v = (4 * o_dual + 2 * c_dual) * d_y;

    SLEQP_CALL(sleqp_vec_push(product, 1, v));
  }

  return SLEQP_OKAY;
}

FuncData* func_data;
SleqpFunc* func;

SleqpParams* params;
SleqpProblem* problem;

void
second_order_setup()
{
  const double inf = sleqp_infinity();

  ASSERT_CALL(sleqp_vec_create(&var_lb, num_variables, num_variables));

  ASSERT_CALL(sleqp_vec_create(&var_ub, num_variables, num_variables));

  for (int i = 0; i < num_variables; ++i)
  {
    ASSERT_CALL(sleqp_vec_push(var_lb, i, -inf));
    ASSERT_CALL(sleqp_vec_push(var_ub, i, inf));
  }

  ASSERT_CALL(sleqp_vec_create(&cons_lb, num_constraints, num_constraints));

  ASSERT_CALL(sleqp_vec_create(&cons_ub, num_constraints, num_constraints));

  ASSERT_CALL(sleqp_vec_push(cons_lb, 0, 0.));
  ASSERT_CALL(sleqp_vec_push(cons_ub, 0, 0.));

  ASSERT_CALL(sleqp_vec_create(&x, num_variables, num_variables));

  ASSERT_CALL(sleqp_vec_push(x, 0, 0.707106781186548));
  ASSERT_CALL(sleqp_vec_push(x, 1, 0.707106781186548));

  ASSERT_CALL(sleqp_malloc(&func_data));

  ASSERT_CALL(sleqp_alloc_array(&func_data->values, num_variables));
  ASSERT_CALL(sleqp_alloc_array(&func_data->duals, num_constraints));
  ASSERT_CALL(sleqp_alloc_array(&func_data->direction, num_variables));

  SleqpFuncCallbacks callbacks = {.set_value = func_set,
                                  .obj_val   = func_obj_val,
                                  .obj_grad  = func_obj_grad,
                                  .cons_val  = func_cons_val,
                                  .cons_jac  = func_cons_jac,
                                  .hess_prod = func_hess_prod,
                                  .func_free = NULL};

  ASSERT_CALL(sleqp_func_create(&func,
                                &callbacks,
                                num_variables,
                                num_constraints,
                                func_data));

  ASSERT_CALL(sleqp_params_create(&params));

  ASSERT_CALL(sleqp_problem_create_simple(&problem,
                                          func,
                                          params,
                                          var_lb,
                                          var_ub,
                                          cons_lb,
                                          cons_ub));
}

void
second_order_teardown()
{
  ASSERT_CALL(sleqp_problem_release(&problem));

  ASSERT_CALL(sleqp_params_release(&params));

  ASSERT_CALL(sleqp_func_release(&func));

  sleqp_free(&func_data->direction);
  sleqp_free(&func_data->duals);
  sleqp_free(&func_data->values);
  sleqp_free(&func_data);

  ASSERT_CALL(sleqp_vec_free(&x));
  ASSERT_CALL(sleqp_vec_free(&cons_ub));
  ASSERT_CALL(sleqp_vec_free(&cons_lb));
  ASSERT_CALL(sleqp_vec_free(&var_ub));
  ASSERT_CALL(sleqp_vec_free(&var_lb));
}

START_TEST(test_second_order_solve)
{
  SleqpVec* expected_solution;

  SleqpParams* params;
  SleqpOptions* options;
  SleqpSolver* solver;

  ASSERT_CALL(sleqp_vec_create(&expected_solution, 2, 2));

  ASSERT_CALL(sleqp_vec_push(expected_solution, 0, 1.));

  ASSERT_CALL(sleqp_params_create(&params));

  ASSERT_CALL(sleqp_options_create(&options));

  ASSERT_CALL(sleqp_solver_create(&solver, problem, params, options, x, NULL));

  // 100 iterations should be plenty...
  ASSERT_CALL(sleqp_solver_solve(solver, 100, -1));

  SleqpIterate* solution_iterate;

  ASSERT_CALL(sleqp_solver_solution(solver, &solution_iterate));

  ck_assert_int_eq(sleqp_solver_status(solver), SLEQP_STATUS_OPTIMAL);

  SleqpVec* actual_solution = sleqp_iterate_primal(solution_iterate);

  ck_assert(sleqp_vec_eq(actual_solution, expected_solution, 1e-6));

  ASSERT_CALL(sleqp_solver_release(&solver));

  ASSERT_CALL(sleqp_options_release(&options));

  ASSERT_CALL(sleqp_params_release(&params));

  ASSERT_CALL(sleqp_vec_free(&expected_solution));
}
END_TEST

Suite*
second_order_test_suite()
{
  Suite* suite;
  TCase* tc_cons;

  suite = suite_create("Constrained tests");

  tc_cons = tcase_create("Constrained solution test");

  tcase_add_checked_fixture(tc_cons, second_order_setup, second_order_teardown);

  tcase_add_test(tc_cons, test_second_order_solve);
  suite_add_tcase(suite, tc_cons);

  return suite;
}

TEST_MAIN(second_order_test_suite)
