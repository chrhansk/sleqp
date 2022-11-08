#include <check.h>
#include <math.h>
#include <stdlib.h>

#include "cmp.h"
#include "lsq.h"
#include "mem.h"
#include "solver.h"

#include "lp/lpi.h"

#include "test_common.h"

#include "rosenbrock_fixture.h"
#include "rosenbrock_lsq_fixture.h"

SleqpOptions* options;
SleqpParams* params;

void
lsq_setup()
{
  rosenbrock_setup();

  rosenbrock_lsq_setup();

  ASSERT_CALL(sleqp_params_create(&params));
  ASSERT_CALL(sleqp_options_create(&options));
}

void
lsq_teardown()
{

  ASSERT_CALL(sleqp_options_release(&options));
  ASSERT_CALL(sleqp_params_release(&params));

  rosenbrock_lsq_teardown();

  rosenbrock_teardown();
}

START_TEST(test_unconstrained_solve)
{
  SleqpProblem* problem;
  SleqpSolver* solver;

  ASSERT_CALL(sleqp_options_set_enum_value(options,
                                           SLEQP_OPTION_ENUM_DERIV_CHECK,
                                           SLEQP_DERIV_CHECK_FIRST));

  ASSERT_CALL(sleqp_problem_create_simple(&problem,
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
                                  rosenbrock_initial,
                                  NULL));

  // 100 iterations should be plenty...
  ASSERT_CALL(sleqp_solver_solve(solver, 100, -1));

  SleqpIterate* solution_iterate;

  ASSERT_CALL(sleqp_solver_solution(solver, &solution_iterate));

  ck_assert_int_eq(sleqp_solver_status(solver), SLEQP_STATUS_OPTIMAL);

  SleqpVec* actual_solution = sleqp_iterate_primal(solution_iterate);

  ck_assert(sleqp_vec_eq(actual_solution, rosenbrock_optimum, 1e-6));

  ASSERT_CALL(sleqp_solver_release(&solver));

  ASSERT_CALL(sleqp_problem_release(&problem));
}
END_TEST

START_TEST(test_lm_factor)
{
  SleqpProblem* problem;
  SleqpSolver* solver;

  ASSERT_CALL(sleqp_params_set_value(params, SLEQP_PARAM_STAT_TOL, 1e-8));

  ASSERT_CALL(sleqp_lsq_func_set_lm_factor(rosenbrock_lsq_func, 1e-2));

  ASSERT_CALL(sleqp_options_set_enum_value(options,
                                           SLEQP_OPTION_ENUM_DERIV_CHECK,
                                           SLEQP_DERIV_CHECK_FIRST));

  ASSERT_CALL(sleqp_problem_create_simple(&problem,
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
                                  rosenbrock_initial,
                                  NULL));

  // 100 iterations should be plenty...
  ASSERT_CALL(sleqp_solver_solve(solver, 100, -1));

  SleqpIterate* solution_iterate;

  ASSERT_CALL(sleqp_solver_solution(solver, &solution_iterate));

  ck_assert_int_eq(sleqp_solver_status(solver), SLEQP_STATUS_OPTIMAL);

  SleqpVec* actual_solution = sleqp_iterate_primal(solution_iterate);

  ck_assert(sleqp_vec_eq(actual_solution, rosenbrock_optimum, 1e-6));

  ASSERT_CALL(sleqp_solver_release(&solver));

  ASSERT_CALL(sleqp_problem_release(&problem));
}
END_TEST

START_TEST(test_lsqr_solve)
{
  SleqpProblem* problem;
  SleqpSolver* solver;

  ASSERT_CALL(sleqp_options_set_enum_value(options,
                                           SLEQP_OPTION_ENUM_TR_SOLVER,
                                           SLEQP_TR_SOLVER_LSQR));

  ASSERT_CALL(sleqp_problem_create_simple(&problem,
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
                                  rosenbrock_initial,
                                  NULL));

  // 100 iterations should be plenty...
  ASSERT_CALL(sleqp_solver_solve(solver, 100, -1));

  SleqpIterate* solution_iterate;

  ASSERT_CALL(sleqp_solver_solution(solver, &solution_iterate));

  ck_assert_int_eq(sleqp_solver_status(solver), SLEQP_STATUS_OPTIMAL);

  SleqpVec* actual_solution = sleqp_iterate_primal(solution_iterate);

  ck_assert(sleqp_vec_eq(actual_solution, rosenbrock_optimum, 1e-6));

  ASSERT_CALL(sleqp_solver_release(&solver));

  ASSERT_CALL(sleqp_problem_release(&problem));
}
END_TEST

START_TEST(test_constrained_lsqr_solve)
{
  SleqpProblem* problem;
  SleqpSolver* solver;

  SleqpMat* linear_coeffs;
  SleqpVec* linear_lb;
  SleqpVec* linear_ub;

  const double inf = sleqp_infinity();

  ASSERT_CALL(sleqp_mat_create(&linear_coeffs, 1, rosenbrock_num_vars, 2));

  ASSERT_CALL(sleqp_mat_push_col(linear_coeffs, 0));

  ASSERT_CALL(sleqp_mat_push(linear_coeffs, 0, 0, 1.));

  ASSERT_CALL(sleqp_mat_push_col(linear_coeffs, 1));

  ASSERT_CALL(sleqp_mat_push(linear_coeffs, 0, 1, 1.));

  ASSERT_CALL(sleqp_vec_create(&linear_lb, 1, 1));
  ASSERT_CALL(sleqp_vec_push(linear_lb, 0, 1.));

  ASSERT_CALL(sleqp_vec_create(&linear_ub, 1, 1));
  ASSERT_CALL(sleqp_vec_push(linear_ub, 0, inf));

  ASSERT_CALL(sleqp_options_set_enum_value(options,
                                           SLEQP_OPTION_ENUM_TR_SOLVER,
                                           SLEQP_TR_SOLVER_LSQR));

  ASSERT_CALL(sleqp_problem_create(&problem,
                                   rosenbrock_lsq_func,
                                   params,
                                   rosenbrock_var_lb,
                                   rosenbrock_var_ub,
                                   rosenbrock_cons_lb,
                                   rosenbrock_cons_ub,
                                   linear_coeffs,
                                   linear_lb,
                                   linear_ub));

  ASSERT_CALL(sleqp_vec_clear(rosenbrock_initial));
  ASSERT_CALL(sleqp_vec_reserve(rosenbrock_initial, 2));

  ASSERT_CALL(sleqp_vec_push(rosenbrock_initial, 0, -1.));
  ASSERT_CALL(sleqp_vec_push(rosenbrock_initial, 1, -1.));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  params,
                                  options,
                                  rosenbrock_initial,
                                  NULL));

  // 100 iterations should be plenty...
  ASSERT_CALL(sleqp_solver_solve(solver, 100, -1));

  SleqpIterate* solution_iterate;

  ASSERT_CALL(sleqp_solver_solution(solver, &solution_iterate));

  ck_assert_int_eq(sleqp_solver_status(solver), SLEQP_STATUS_OPTIMAL);

  SleqpVec* actual_solution = sleqp_iterate_primal(solution_iterate);

  ck_assert(sleqp_vec_eq(actual_solution, rosenbrock_optimum, 1e-6));

  ASSERT_CALL(sleqp_solver_release(&solver));

  ASSERT_CALL(sleqp_problem_release(&problem));

  ASSERT_CALL(sleqp_vec_free(&linear_ub));
  ASSERT_CALL(sleqp_vec_free(&linear_lb));

  ASSERT_CALL(sleqp_mat_release(&linear_coeffs));
}
END_TEST

START_TEST(test_scaled_solve)
{
  SleqpProblem* problem;
  SleqpScaling* scaling;
  SleqpSolver* solver;

  ASSERT_CALL(
    sleqp_scaling_create(&scaling, rosenbrock_num_vars, rosenbrock_num_cons));

  ASSERT_CALL(sleqp_scaling_set_obj_weight(scaling, -2));

  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling, 0, -1));

  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling, 1, -1));

  ASSERT_CALL(sleqp_options_set_enum_value(options,
                                           SLEQP_OPTION_ENUM_DERIV_CHECK,
                                           SLEQP_DERIV_CHECK_FIRST));

  ASSERT_CALL(sleqp_problem_create_simple(&problem,
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
                                  rosenbrock_initial,
                                  scaling));

  // 100 iterations should be plenty...
  ASSERT_CALL(sleqp_solver_solve(solver, 100, -1));

  SleqpIterate* solution_iterate;

  ASSERT_CALL(sleqp_solver_solution(solver, &solution_iterate));

  ck_assert_int_eq(sleqp_solver_status(solver), SLEQP_STATUS_OPTIMAL);

  SleqpVec* actual_solution = sleqp_iterate_primal(solution_iterate);

  ck_assert(sleqp_vec_eq(actual_solution, rosenbrock_optimum, 1e-6));

  ASSERT_CALL(sleqp_solver_release(&solver));

  ASSERT_CALL(sleqp_problem_release(&problem));

  ASSERT_CALL(sleqp_scaling_release(&scaling));
}
END_TEST

Suite*
lsq_test_suite()
{
  Suite* suite;
  TCase* tc_uncons;

  suite = suite_create("LSQ tests");

  tc_uncons = tcase_create("LSQ solution test");

  tcase_add_checked_fixture(tc_uncons, lsq_setup, lsq_teardown);

  tcase_add_test(tc_uncons, test_unconstrained_solve);
  tcase_add_test(tc_uncons, test_lm_factor);
  tcase_add_test(tc_uncons, test_lsqr_solve);
  tcase_add_test(tc_uncons, test_constrained_lsqr_solve);
  tcase_add_test(tc_uncons, test_scaled_solve);

  suite_add_tcase(suite, tc_uncons);

  return suite;
}

TEST_MAIN(lsq_test_suite)
