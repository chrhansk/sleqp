#include <check.h>
#include <stdlib.h>

#include "cmp.h"
#include "direction.h"
#include "gauss_newton.h"
#include "problem.h"
#include "sparse/pub_vec.h"
#include "util.h"
#include "working_step.h"

#include "aug_jac/unconstrained_aug_jac.h"

#include "test_common.h"
#include "zero_func.h"

#include "linear_lsq.h"

SleqpProblem* problem;

SleqpWorkingStep* working_step;

SleqpEQPSolver* gauss_newton_solver;

SleqpIterate* iterate;

SleqpAugJac* aug_jac;

SleqpDirection* direction;
SleqpVec* cons_dual;
SleqpVec* point;
SleqpVec* initial;

SleqpFunc* zero_func;

SleqpVec* zero_vec;

void
unconstrained_setup()
{
  linear_lsq_setup();

  ASSERT_CALL(sleqp_problem_create_simple(&problem,
                                          linear_lsq_func,
                                          linear_lsq_settings,
                                          linear_lsq_var_lb,
                                          linear_lsq_var_ub,
                                          linear_lsq_cons_lb,
                                          linear_lsq_cons_ub));

  ASSERT_CALL(
    sleqp_working_step_create(&working_step, problem, linear_lsq_settings));

  ASSERT_CALL(sleqp_gauss_newton_solver_create(&gauss_newton_solver,
                                               problem,
                                               linear_lsq_settings,
                                               working_step));

  ASSERT_CALL(sleqp_iterate_create(&iterate, problem, linear_lsq_initial));

  ASSERT_CALL(sleqp_unconstrained_aug_jac_create(&aug_jac, problem));

  ASSERT_CALL(sleqp_direction_create(&direction, problem, linear_lsq_settings));

  ASSERT_CALL(sleqp_vec_create_empty(&cons_dual, linear_lsq_num_constraints));

  ASSERT_CALL(sleqp_vec_create_empty(&point, linear_lsq_num_variables));

  ASSERT_CALL(sleqp_vec_create_full(&initial, linear_lsq_num_variables));
}

void
unconstrained_teardown()
{
  ASSERT_CALL(sleqp_vec_free(&initial));

  ASSERT_CALL(sleqp_vec_free(&point));

  ASSERT_CALL(sleqp_vec_free(&cons_dual));

  ASSERT_CALL(sleqp_direction_release(&direction));

  ASSERT_CALL(sleqp_aug_jac_release(&aug_jac));

  ASSERT_CALL(sleqp_iterate_release(&iterate));

  ASSERT_CALL(sleqp_eqp_solver_release(&gauss_newton_solver));

  ASSERT_CALL(sleqp_working_step_release(&working_step));

  ASSERT_CALL(sleqp_problem_release(&problem));

  linear_lsq_teardown();
}

void
compute_point(const SleqpVec* initial, SleqpVec* point)
{
  const double penalty_parameter = 1.;
  const double trust_radius      = 1e6;

  const double zero_eps
    = sleqp_settings_real_value(linear_lsq_settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  ASSERT_CALL(sleqp_vec_copy(initial, sleqp_iterate_primal(iterate)));

  bool reject;

  ASSERT_CALL(
    sleqp_set_and_evaluate(problem, iterate, SLEQP_VALUE_REASON_NONE, &reject));

  assert(!reject);

  ASSERT_CALL(sleqp_eqp_solver_set_iterate(gauss_newton_solver,
                                           iterate,
                                           aug_jac,
                                           trust_radius,
                                           penalty_parameter));

  ASSERT_CALL(sleqp_eqp_solver_compute_direction(gauss_newton_solver,
                                                 cons_dual,
                                                 direction));

  SleqpVec* step = sleqp_direction_primal(direction);

  ASSERT_CALL(sleqp_vec_add(initial, step, zero_eps, point));
}

// LSQR direction should point to optimum
START_TEST(test_unconstrained_solve)
{
  const double eps = sleqp_settings_real_value(linear_lsq_settings, SLEQP_SETTINGS_REAL_EPS);

  compute_point(linear_lsq_initial, point);

  ck_assert(sleqp_vec_eq(point, linear_lsq_optimal, eps));
}
END_TEST

// LSQR direction should point to optimum
START_TEST(test_unconstrained_start_high)
{
  const double eps = sleqp_settings_real_value(linear_lsq_settings, SLEQP_SETTINGS_REAL_EPS);

  {
    double values[] = {20., 20.};

    ASSERT_CALL(
      sleqp_vec_set_from_raw(initial, values, linear_lsq_num_variables, 0.));
  }

  compute_point(initial, point);

  ck_assert(sleqp_vec_eq(point, linear_lsq_optimal, eps));
}
END_TEST

// LSQR direction should point to optimum
START_TEST(test_unconstrained_start_low)
{
  const double eps = sleqp_settings_real_value(linear_lsq_settings, SLEQP_SETTINGS_REAL_EPS);

  {
    double values[] = {-10., -10.};

    ASSERT_CALL(
      sleqp_vec_set_from_raw(initial, values, linear_lsq_num_variables, 0.));
  }

  compute_point(initial, point);

  ck_assert(sleqp_vec_eq(point, linear_lsq_optimal, eps));
}
END_TEST

// Insufficient trust radius, solution should be on the boundary
START_TEST(test_unconstrained_small)
{
  const double penalty_parameter = 1.;
  const double trust_radius      = 1.;

  const double eps = sleqp_settings_real_value(linear_lsq_settings, SLEQP_SETTINGS_REAL_EPS);

  ASSERT_CALL(sleqp_eqp_solver_set_iterate(gauss_newton_solver,
                                           iterate,
                                           aug_jac,
                                           trust_radius,
                                           penalty_parameter));

  ASSERT_CALL(sleqp_eqp_solver_compute_direction(gauss_newton_solver,
                                                 cons_dual,
                                                 direction));

  SleqpVec* step = sleqp_direction_primal(direction);

  ck_assert(sleqp_is_eq(sleqp_vec_norm(step), trust_radius, eps));
}
END_TEST

void
constrained_setup()
{
  linear_lsq_setup();

  int num_constraints = 0;
  int num_residuals   = 0;

  ASSERT_CALL(zero_lsq_func_create(&zero_func,
                                   linear_lsq_settings,
                                   linear_lsq_num_variables,
                                   num_constraints,
                                   num_residuals));

  ASSERT_CALL(sleqp_vec_create_empty(&zero_vec, 0));

  ASSERT_CALL(sleqp_problem_create(&problem,
                                   zero_func,
                                   linear_lsq_settings,
                                   linear_lsq_var_lb,
                                   linear_lsq_var_ub,
                                   zero_vec,
                                   zero_vec,
                                   linear_lsq_matrix,
                                   linear_lsq_rhs,
                                   linear_lsq_rhs));

  ASSERT_CALL(
    sleqp_working_step_create(&working_step, problem, linear_lsq_settings));

  ASSERT_CALL(sleqp_gauss_newton_solver_create(&gauss_newton_solver,
                                               problem,
                                               linear_lsq_settings,
                                               working_step));

  ASSERT_CALL(sleqp_iterate_create(&iterate, problem, linear_lsq_initial));

  ASSERT_CALL(sleqp_unconstrained_aug_jac_create(&aug_jac, problem));

  ASSERT_CALL(sleqp_direction_create(&direction, problem, linear_lsq_settings));

  ASSERT_CALL(sleqp_vec_create_empty(&cons_dual, linear_lsq_num_constraints));

  ASSERT_CALL(sleqp_vec_create_empty(&point, linear_lsq_num_variables));

  ASSERT_CALL(sleqp_vec_create_full(&initial, linear_lsq_num_variables));
}

void
constrained_teardown()
{
  ASSERT_CALL(sleqp_vec_free(&initial));

  ASSERT_CALL(sleqp_vec_free(&point));

  ASSERT_CALL(sleqp_vec_free(&cons_dual));

  ASSERT_CALL(sleqp_direction_release(&direction));

  ASSERT_CALL(sleqp_aug_jac_release(&aug_jac));

  ASSERT_CALL(sleqp_iterate_release(&iterate));

  ASSERT_CALL(sleqp_eqp_solver_release(&gauss_newton_solver));

  ASSERT_CALL(sleqp_working_step_release(&working_step));

  ASSERT_CALL(sleqp_problem_release(&problem));

  ASSERT_CALL(sleqp_vec_free(&zero_vec));

  ASSERT_CALL(sleqp_func_release(&zero_func));

  linear_lsq_teardown();
}

START_TEST(test_constrained_start_high)
{
  const double eps = sleqp_settings_real_value(linear_lsq_settings, SLEQP_SETTINGS_REAL_EPS);

  {
    double values[] = {20., 20.};

    ASSERT_CALL(
      sleqp_vec_set_from_raw(initial, values, linear_lsq_num_variables, 0.));
  }

  compute_point(initial, point);

  ck_assert(sleqp_vec_eq(point, linear_lsq_optimal, eps));
}
END_TEST

START_TEST(test_constrained_start_low)
{
  const double eps = sleqp_settings_real_value(linear_lsq_settings, SLEQP_SETTINGS_REAL_EPS);

  {
    double values[] = {-10., -10.};

    ASSERT_CALL(
      sleqp_vec_set_from_raw(initial, values, linear_lsq_num_variables, 0.));
  }

  compute_point(initial, point);

  ck_assert(sleqp_vec_eq(point, linear_lsq_optimal, eps));
}
END_TEST

Suite*
gauss_newton_test_suite()
{
  Suite* suite;
  TCase* tc_solve_uncons;
  TCase* tc_solve_cons;

  suite = suite_create("Gauss-Newton tests");

  tc_solve_uncons = tcase_create("Unconstrained solutions");

  tcase_add_checked_fixture(tc_solve_uncons,
                            unconstrained_setup,
                            unconstrained_teardown);

  tcase_add_test(tc_solve_uncons, test_unconstrained_solve);

  tcase_add_test(tc_solve_uncons, test_unconstrained_start_high);

  tcase_add_test(tc_solve_uncons, test_unconstrained_start_low);

  tcase_add_test(tc_solve_uncons, test_unconstrained_small);

  suite_add_tcase(suite, tc_solve_uncons);

  tc_solve_cons = tcase_create("Constrained solutions");

  tcase_add_checked_fixture(tc_solve_cons,
                            constrained_setup,
                            constrained_teardown);

  tcase_add_test(tc_solve_cons, test_constrained_start_high);

  tcase_add_test(tc_solve_cons, test_constrained_start_low);

  suite_add_tcase(suite, tc_solve_cons);

  return suite;
}

TEST_MAIN(gauss_newton_test_suite)
