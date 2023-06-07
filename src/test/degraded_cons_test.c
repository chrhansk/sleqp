#include <check.h>
#include <stdlib.h>

#include "test_common.h"

#include "sleqp.h"

/*
 * Instance with two linearly dependent (linear) constraints,
 * not satisfying LICQ / MFCQ constraint qualifications
 */

SleqpSettings* settings;
SleqpProblem* problem;

#define DEGRADED_NUM_VARS 2
#define DEGRADED_NUM_CONS 0

#define DEGRADED_NUM_LINEAR 3

double values[DEGRADED_NUM_VARS];
const double target[DEGRADED_NUM_VARS] = {1., 2.5};

double
sq(double x)
{
  return x * x;
}

SLEQP_RETCODE
degraded_set(SleqpFunc* func,
             SleqpVec* x,
             SLEQP_VALUE_REASON reason,
             bool* reject,
             void* func_data)
{
  SLEQP_CALL(sleqp_vec_to_raw(x, values));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
degraded_obj_val(SleqpFunc* func, double* obj_val, void* func_data)
{
  *obj_val = .5 * (sq(values[0] - target[0]) + sq(values[1] - target[1]));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
degraded_obj_grad(SleqpFunc* func, SleqpVec* obj_grad, void* func_data)
{
  SLEQP_CALL(sleqp_vec_push(obj_grad, 0, values[0] - target[0]));
  SLEQP_CALL(sleqp_vec_push(obj_grad, 1, values[1] - target[1]));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
degraded_hess_prod(SleqpFunc* func,
                   const SleqpVec* direction,
                   const SleqpVec* cons_duals,
                   SleqpVec* product,
                   void* func_data)
{
  SLEQP_CALL(sleqp_vec_reserve(product, DEGRADED_NUM_VARS));

  SLEQP_CALL(sleqp_vec_copy(direction, product));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
degraded_func_create(SleqpFunc** fstar)
{
  SleqpFuncCallbacks callbacks = {.set_value = degraded_set,
                                  .obj_val   = degraded_obj_val,
                                  .obj_grad  = degraded_obj_grad,
                                  .cons_val  = NULL,
                                  .cons_jac  = NULL,
                                  .hess_prod = degraded_hess_prod,
                                  .func_free = NULL};

  SLEQP_CALL(sleqp_func_create(fstar,
                               &callbacks,
                               DEGRADED_NUM_VARS,
                               DEGRADED_NUM_CONS,
                               NULL));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
linear_cons_create(SleqpMat** cstar,
                   SleqpVec** lbstar,
                   SleqpVec** ubstar,
                   bool swapped)
{
  const double inf = sleqp_infinity();

  SLEQP_CALL(sleqp_vec_create_full(lbstar, DEGRADED_NUM_LINEAR));
  SLEQP_CALL(sleqp_vec_create_full(ubstar, DEGRADED_NUM_LINEAR));

  const int nnz = DEGRADED_NUM_LINEAR * DEGRADED_NUM_VARS;

  SLEQP_CALL(
    sleqp_mat_create(cstar, DEGRADED_NUM_LINEAR, DEGRADED_NUM_VARS, nnz));

  SleqpVec* linear_lb     = *lbstar;
  SleqpVec* linear_ub     = *ubstar;
  SleqpMat* linear_coeffs = *cstar;

  if (swapped)
  {
    SLEQP_CALL(sleqp_vec_push(linear_lb, 0, 2.));
    SLEQP_CALL(sleqp_vec_push(linear_lb, 1, -6.));
    SLEQP_CALL(sleqp_vec_push(linear_lb, 2, -2.));
  }
  else
  {
    SLEQP_CALL(sleqp_vec_push(linear_lb, 0, -2.));
    SLEQP_CALL(sleqp_vec_push(linear_lb, 1, -6.));
    SLEQP_CALL(sleqp_vec_push(linear_lb, 2, 2.));
  }

  SLEQP_CALL(sleqp_vec_fill(linear_ub, inf));

  if (swapped)
  {
    SLEQP_CALL(sleqp_mat_push(linear_coeffs, 0, 0, -1.));
    SLEQP_CALL(sleqp_mat_push(linear_coeffs, 1, 0, -1.));
    SLEQP_CALL(sleqp_mat_push(linear_coeffs, 2, 0, 1.));

    SLEQP_CALL(sleqp_mat_push_col(linear_coeffs, 1));

    SLEQP_CALL(sleqp_mat_push(linear_coeffs, 0, 1, 2.));
    SLEQP_CALL(sleqp_mat_push(linear_coeffs, 1, 1, -2.));
    SLEQP_CALL(sleqp_mat_push(linear_coeffs, 2, 1, -2.));
  }
  else
  {
    SLEQP_CALL(sleqp_mat_push(linear_coeffs, 0, 0, 1.));
    SLEQP_CALL(sleqp_mat_push(linear_coeffs, 1, 0, -1.));
    SLEQP_CALL(sleqp_mat_push(linear_coeffs, 2, 0, -1.));

    SLEQP_CALL(sleqp_mat_push_col(linear_coeffs, 1));

    SLEQP_CALL(sleqp_mat_push(linear_coeffs, 0, 1, -2.));
    SLEQP_CALL(sleqp_mat_push(linear_coeffs, 1, 1, -2.));
    SLEQP_CALL(sleqp_mat_push(linear_coeffs, 2, 1, 2.));
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
degraded_problem_create(SleqpProblem** star, SleqpSettings* settings, bool swapped)
{
  SleqpFunc* func;

  SleqpVec* var_lb;
  SleqpVec* var_ub;

  SleqpVec* cons_bounds;

  SleqpVec* linear_lb;
  SleqpVec* linear_ub;
  SleqpMat* linear_coeffs;

  const double inf = sleqp_infinity();

  SLEQP_CALL(degraded_func_create(&func));

  SLEQP_CALL(sleqp_vec_create_full(&var_lb, DEGRADED_NUM_VARS));
  SLEQP_CALL(sleqp_vec_create_full(&var_ub, DEGRADED_NUM_VARS));

  SLEQP_CALL(sleqp_vec_fill(var_lb, -inf));
  SLEQP_CALL(sleqp_vec_fill(var_ub, inf));

  SLEQP_CALL(sleqp_vec_create_empty(&cons_bounds, 0.));

  SLEQP_CALL(
    linear_cons_create(&linear_coeffs, &linear_lb, &linear_ub, swapped));

  SLEQP_CALL(sleqp_problem_create(star,
                                  func,
                                  var_lb,
                                  var_ub,
                                  cons_bounds,
                                  cons_bounds,
                                  linear_coeffs,
                                  linear_lb,
                                  linear_ub,
                                  settings));

  SLEQP_CALL(sleqp_mat_release(&linear_coeffs));
  SLEQP_CALL(sleqp_vec_free(&linear_ub));
  SLEQP_CALL(sleqp_vec_free(&linear_lb));

  SLEQP_CALL(sleqp_vec_free(&cons_bounds));

  SLEQP_CALL(sleqp_vec_free(&var_ub));
  SLEQP_CALL(sleqp_vec_free(&var_lb));

  SLEQP_CALL(sleqp_func_release(&func));

  return SLEQP_OKAY;
}

SleqpSettings* settings;
SleqpVec* initial;

void
degraded_setup()
{
  ASSERT_CALL(sleqp_settings_create(&settings));

  ASSERT_CALL(sleqp_settings_set_enum_value(
    settings,
    SLEQP_SETTINGS_ENUM_DERIV_CHECK,
    SLEQP_DERIV_CHECK_FIRST | SLEQP_DERIV_CHECK_SECOND_EXHAUSTIVE));

  ASSERT_CALL(sleqp_vec_create_full(&initial, DEGRADED_NUM_VARS));

  ASSERT_CALL(sleqp_vec_push(initial, 0, 0.));
  ASSERT_CALL(sleqp_vec_push(initial, 1, 0.));
}

void
degraded_teardown()
{
  ASSERT_CALL(sleqp_vec_free(&initial));

  ASSERT_CALL(sleqp_settings_release(&settings));
}

START_TEST(test_solve)
{
  SleqpProblem* problem;
  SleqpSolver* solver;

  ASSERT_CALL(degraded_problem_create(&problem, settings, false));

  ASSERT_CALL(
    sleqp_solver_create(&solver, problem, initial, NULL));

  ASSERT_CALL(sleqp_solver_solve(solver, 100, SLEQP_NONE));

  SleqpIterate* iterate;

  ASSERT_CALL(sleqp_solver_solution(solver, &iterate));

  ck_assert_int_eq(sleqp_solver_status(solver), SLEQP_STATUS_OPTIMAL);

  ASSERT_CALL(sleqp_solver_release(&solver));

  ASSERT_CALL(sleqp_problem_release(&problem));
}
END_TEST

START_TEST(test_solve_swapped)
{
  SleqpProblem* problem;
  SleqpSolver* solver;

  ASSERT_CALL(degraded_problem_create(&problem, settings, true));

  ASSERT_CALL(
    sleqp_solver_create(&solver, problem, initial, NULL));

  ASSERT_CALL(sleqp_solver_solve(solver, 10, SLEQP_NONE));

  ck_assert_int_eq(sleqp_solver_status(solver), SLEQP_STATUS_OPTIMAL);

  ASSERT_CALL(sleqp_solver_release(&solver));

  ASSERT_CALL(sleqp_problem_release(&problem));
}
END_TEST

Suite*
degraded_test_suite()
{
  Suite* suite;
  TCase* tc_deg;

  suite = suite_create("Degraded tests");

  tc_deg = tcase_create("Degraded solution test");

  tcase_add_checked_fixture(tc_deg, degraded_setup, degraded_teardown);

  tcase_add_test(tc_deg, test_solve);

  tcase_add_test(tc_deg, test_solve_swapped);

  suite_add_tcase(suite, tc_deg);

  return suite;
}

TEST_MAIN(degraded_test_suite)
