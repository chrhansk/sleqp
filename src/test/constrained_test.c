#include <check.h>
#include <stdlib.h>

#include "cmp.h"
#include "mem.h"
#include "solver.h"
#include "util.h"

#include "constrained_fixture.h"
#include "test_common.h"

SleqpParams* params;
SleqpOptions* options;
SleqpProblem* problem;

void
constrained_test_setup()
{
  constrained_setup();

  ASSERT_CALL(sleqp_params_create(&params));

  ASSERT_CALL(sleqp_options_create(&options));

  ASSERT_CALL(sleqp_problem_create_simple(&problem,
                                          constrained_func,
                                          params,
                                          constrained_var_lb,
                                          constrained_var_ub,
                                          constrained_cons_lb,
                                          constrained_cons_ub));
}

void
constrained_test_teardown()
{
  ASSERT_CALL(sleqp_problem_release(&problem));

  ASSERT_CALL(sleqp_options_release(&options));

  ASSERT_CALL(sleqp_params_release(&params));

  constrained_teardown();
}

START_TEST(test_solve)
{
  SleqpSolver* solver;

  ASSERT_CALL(sleqp_options_set_int(options,
                                    SLEQP_OPTION_INT_DERIV_CHECK,
                                    SLEQP_DERIV_CHECK_FIRST
                                      | SLEQP_DERIV_CHECK_SECOND_EXHAUSTIVE));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  params,
                                  options,
                                  constrained_initial,
                                  NULL));

  // 100 iterations should be plenty...
  ASSERT_CALL(sleqp_solver_solve(solver, 100, -1));

  SleqpIterate* solution_iterate;

  ASSERT_CALL(sleqp_solver_get_solution(solver, &solution_iterate));

  ck_assert_int_eq(sleqp_solver_get_status(solver), SLEQP_STATUS_OPTIMAL);

  SleqpSparseVec* actual_solution = sleqp_iterate_get_primal(solution_iterate);

  ck_assert(sleqp_sparse_vector_eq(actual_solution, constrained_optimal, 1e-6));

  ASSERT_CALL(sleqp_solver_release(&solver));
}
END_TEST

START_TEST(test_exact_linesearch)
{
  SleqpSolver* solver;

  ASSERT_CALL(sleqp_options_set_int(options,
                                    SLEQP_OPTION_INT_LINESEARCH,
                                    SLEQP_LINESEARCH_EXACT));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  params,
                                  options,
                                  constrained_initial,
                                  NULL));

  // 100 iterations should be plenty...
  ASSERT_CALL(sleqp_solver_solve(solver, 100, -1));

  SleqpIterate* solution_iterate;

  ASSERT_CALL(sleqp_solver_get_solution(solver, &solution_iterate));

  ck_assert_int_eq(sleqp_solver_get_status(solver), SLEQP_STATUS_OPTIMAL);

  SleqpSparseVec* actual_solution = sleqp_iterate_get_primal(solution_iterate);

  ck_assert(sleqp_sparse_vector_eq(actual_solution, constrained_optimal, 1e-6));

  ASSERT_CALL(sleqp_solver_release(&solver));
}
END_TEST

START_TEST(test_parametric_solve)
{
  SleqpSolver* solver;

  ASSERT_CALL(sleqp_options_set_int(options,
                                    SLEQP_OPTION_INT_PARAMETRIC_CAUCHY,
                                    SLEQP_PARAMETRIC_CAUCHY_COARSE));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  params,
                                  options,
                                  constrained_initial,
                                  NULL));

  // 100 iterations should be plenty...
  ASSERT_CALL(sleqp_solver_solve(solver, 100, -1));

  SleqpIterate* solution_iterate;

  ASSERT_CALL(sleqp_solver_get_solution(solver, &solution_iterate));

  ck_assert_int_eq(sleqp_solver_get_status(solver), SLEQP_STATUS_OPTIMAL);

  SleqpSparseVec* actual_solution = sleqp_iterate_get_primal(solution_iterate);

  ck_assert(sleqp_sparse_vector_eq(actual_solution, constrained_optimal, 1e-6));

  ASSERT_CALL(sleqp_solver_release(&solver));
}
END_TEST

START_TEST(test_sr1_solve)
{
  ASSERT_CALL(sleqp_options_set_int(options,
                                    SLEQP_OPTION_INT_DERIV_CHECK,
                                    SLEQP_DERIV_CHECK_FIRST));

  ASSERT_CALL(sleqp_options_set_int(options,
                                    SLEQP_OPTION_INT_HESSIAN_EVAL,
                                    SLEQP_HESSIAN_EVAL_SR1));

  SleqpSolver* solver;

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  params,
                                  options,
                                  constrained_initial,
                                  NULL));

  // 100 iterations should be plenty...
  ASSERT_CALL(sleqp_solver_solve(solver, 100, -1));

  SleqpIterate* solution_iterate;

  ASSERT_CALL(sleqp_solver_get_solution(solver, &solution_iterate));

  ck_assert_int_eq(sleqp_solver_get_status(solver), SLEQP_STATUS_OPTIMAL);

  SleqpSparseVec* actual_solution = sleqp_iterate_get_primal(solution_iterate);

  ck_assert(sleqp_sparse_vector_eq(actual_solution, constrained_optimal, 1e-6));

  ASSERT_CALL(sleqp_solver_release(&solver));
}
END_TEST

START_TEST(test_bfgs_solve_no_sizing)
{
  ASSERT_CALL(sleqp_options_set_int(options,
                                    SLEQP_OPTION_INT_DERIV_CHECK,
                                    SLEQP_DERIV_CHECK_FIRST));

  ASSERT_CALL(sleqp_options_set_int(options,
                                    SLEQP_OPTION_INT_BFGS_SIZING,
                                    SLEQP_BFGS_SIZING_NONE));

  ASSERT_CALL(sleqp_options_set_int(options,
                                    SLEQP_OPTION_INT_HESSIAN_EVAL,
                                    SLEQP_HESSIAN_EVAL_DAMPED_BFGS));

  SleqpSolver* solver;

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  params,
                                  options,
                                  constrained_initial,
                                  NULL));

  // 100 iterations should be plenty...
  ASSERT_CALL(sleqp_solver_solve(solver, 100, -1));

  SleqpIterate* solution_iterate;

  ASSERT_CALL(sleqp_solver_get_solution(solver, &solution_iterate));

  ck_assert_int_eq(sleqp_solver_get_status(solver), SLEQP_STATUS_OPTIMAL);

  SleqpSparseVec* actual_solution = sleqp_iterate_get_primal(solution_iterate);

  ck_assert(sleqp_sparse_vector_eq(actual_solution, constrained_optimal, 1e-6));

  ASSERT_CALL(sleqp_solver_release(&solver));
}
END_TEST

START_TEST(test_bfgs_solve_centered_ol_sizing)
{
  ASSERT_CALL(sleqp_options_set_int(options,
                                    SLEQP_OPTION_INT_BFGS_SIZING,
                                    SLEQP_BFGS_SIZING_CENTERED_OL));

  ASSERT_CALL(sleqp_options_set_int(options,
                                    SLEQP_OPTION_INT_HESSIAN_EVAL,
                                    SLEQP_HESSIAN_EVAL_DAMPED_BFGS));

  SleqpSolver* solver;

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  params,
                                  options,
                                  constrained_initial,
                                  NULL));

  // 100 iterations should be plenty...
  ASSERT_CALL(sleqp_solver_solve(solver, 100, -1));

  SleqpIterate* solution_iterate;

  ASSERT_CALL(sleqp_solver_get_solution(solver, &solution_iterate));

  ck_assert_int_eq(sleqp_solver_get_status(solver), SLEQP_STATUS_OPTIMAL);

  SleqpSparseVec* actual_solution = sleqp_iterate_get_primal(solution_iterate);

  ck_assert(sleqp_sparse_vector_eq(actual_solution, constrained_optimal, 1e-6));

  ASSERT_CALL(sleqp_solver_release(&solver));
}
END_TEST

START_TEST(test_unscaled_solve)
{
  ASSERT_CALL(sleqp_options_set_int(options,
                                    SLEQP_OPTION_INT_DERIV_CHECK,
                                    SLEQP_DERIV_CHECK_FIRST));

  SleqpSolver* solver;

  SleqpScaling* scaling;

  ASSERT_CALL(sleqp_scaling_create(&scaling,
                                   constrained_num_variables,
                                   constrained_num_constraints));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  params,
                                  options,
                                  constrained_initial,
                                  scaling));

  // 100 iterations should be plenty...
  ASSERT_CALL(sleqp_solver_solve(solver, 100, -1));

  SleqpIterate* solution_iterate;

  ASSERT_CALL(sleqp_solver_get_solution(solver, &solution_iterate));

  ck_assert_int_eq(sleqp_solver_get_status(solver), SLEQP_STATUS_OPTIMAL);

  SleqpSparseVec* actual_solution = sleqp_iterate_get_primal(solution_iterate);

  ck_assert(sleqp_sparse_vector_eq(actual_solution, constrained_optimal, 1e-6));

  ASSERT_CALL(sleqp_solver_release(&solver));

  ASSERT_CALL(sleqp_scaling_release(&scaling));
}
END_TEST

START_TEST(test_scaled_solve)
{
  ASSERT_CALL(sleqp_options_set_int(options,
                                    SLEQP_OPTION_INT_DERIV_CHECK,
                                    SLEQP_DERIV_CHECK_FIRST
                                      | SLEQP_DERIV_CHECK_SECOND_EXHAUSTIVE));

  SleqpSolver* solver;

  SleqpScaling* scaling;

  ASSERT_CALL(sleqp_scaling_create(&scaling,
                                   constrained_num_variables,
                                   constrained_num_constraints));

  ASSERT_CALL(sleqp_scaling_set_func_weight(scaling, 2));

  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling, 0, -5));
  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling, 1, 5));

  ASSERT_CALL(sleqp_scaling_set_cons_weight(scaling, 0, -1));
  ASSERT_CALL(sleqp_scaling_set_cons_weight(scaling, 1, 2));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  params,
                                  options,
                                  constrained_initial,
                                  scaling));

  // 100 iterations should be plenty...
  ASSERT_CALL(sleqp_solver_solve(solver, 1000, -1));

  SleqpIterate* solution_iterate;

  ASSERT_CALL(sleqp_solver_get_solution(solver, &solution_iterate));

  ck_assert_int_eq(sleqp_solver_get_status(solver), SLEQP_STATUS_OPTIMAL);

  SleqpSparseVec* actual_solution = sleqp_iterate_get_primal(solution_iterate);

  ck_assert(sleqp_sparse_vector_eq(actual_solution, constrained_optimal, 1e-6));

  ASSERT_CALL(sleqp_solver_release(&solver));

  ASSERT_CALL(sleqp_scaling_release(&scaling));
}
END_TEST

START_TEST(test_scaled_sr1_solve)
{
  ASSERT_CALL(sleqp_options_set_int(options,
                                    SLEQP_OPTION_INT_DERIV_CHECK,
                                    SLEQP_DERIV_CHECK_FIRST));

  SleqpSolver* solver;

  SleqpScaling* scaling;

  ASSERT_CALL(sleqp_options_set_int(options,
                                    SLEQP_OPTION_INT_HESSIAN_EVAL,
                                    SLEQP_HESSIAN_EVAL_SR1));

  ASSERT_CALL(sleqp_scaling_create(&scaling,
                                   constrained_num_variables,
                                   constrained_num_constraints));

  ASSERT_CALL(sleqp_scaling_set_func_weight(scaling, 2));

  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling, 0, -5));
  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling, 1, 5));

  ASSERT_CALL(sleqp_scaling_set_cons_weight(scaling, 0, -1));
  ASSERT_CALL(sleqp_scaling_set_cons_weight(scaling, 1, 2));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  params,
                                  options,
                                  constrained_initial,
                                  scaling));

  // 100 iterations should be plenty...
  ASSERT_CALL(sleqp_solver_solve(solver, 1000, -1));

  SleqpIterate* solution_iterate;

  ASSERT_CALL(sleqp_solver_get_solution(solver, &solution_iterate));

  ck_assert_int_eq(sleqp_solver_get_status(solver), SLEQP_STATUS_OPTIMAL);

  SleqpSparseVec* actual_solution = sleqp_iterate_get_primal(solution_iterate);

  ck_assert(sleqp_sparse_vector_eq(actual_solution, constrained_optimal, 1e-6));

  ASSERT_CALL(sleqp_solver_release(&solver));

  ASSERT_CALL(sleqp_scaling_release(&scaling));
}
END_TEST

START_TEST(test_scaled_bfgs_solve)
{
  ASSERT_CALL(sleqp_options_set_int(options,
                                    SLEQP_OPTION_INT_DERIV_CHECK,
                                    SLEQP_DERIV_CHECK_FIRST));

  SleqpSolver* solver;

  SleqpScaling* scaling;

  ASSERT_CALL(sleqp_options_set_int(options,
                                    SLEQP_OPTION_INT_HESSIAN_EVAL,
                                    SLEQP_HESSIAN_EVAL_DAMPED_BFGS));

  ASSERT_CALL(sleqp_scaling_create(&scaling,
                                   constrained_num_variables,
                                   constrained_num_constraints));

  ASSERT_CALL(sleqp_scaling_set_func_weight(scaling, 2));

  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling, 0, -5));
  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling, 1, 5));

  ASSERT_CALL(sleqp_scaling_set_cons_weight(scaling, 0, -1));
  ASSERT_CALL(sleqp_scaling_set_cons_weight(scaling, 1, 2));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  params,
                                  options,
                                  constrained_initial,
                                  scaling));

  // 100 iterations should be plenty...
  ASSERT_CALL(sleqp_solver_solve(solver, 1000, -1));

  SleqpIterate* solution_iterate;

  ASSERT_CALL(sleqp_solver_get_solution(solver, &solution_iterate));

  ck_assert_int_eq(sleqp_solver_get_status(solver), SLEQP_STATUS_OPTIMAL);

  SleqpSparseVec* actual_solution = sleqp_iterate_get_primal(solution_iterate);

  ck_assert(sleqp_sparse_vector_eq(actual_solution, constrained_optimal, 1e-6));

  ASSERT_CALL(sleqp_solver_release(&solver));

  ASSERT_CALL(sleqp_scaling_release(&scaling));
}
END_TEST

START_TEST(test_auto_scaled_solve)
{
  ASSERT_CALL(sleqp_options_set_int(options,
                                    SLEQP_OPTION_INT_DERIV_CHECK,
                                    SLEQP_DERIV_CHECK_FIRST));

  SleqpSolver* solver;

  SleqpScaling* scaling;

  SleqpIterate* iterate;

  const double eps = sleqp_params_get(params, SLEQP_PARAM_EPS);

  ASSERT_CALL(sleqp_iterate_create(&iterate, problem, constrained_initial));

  ASSERT_CALL(
    sleqp_set_and_evaluate(problem, iterate, SLEQP_VALUE_REASON_NONE, NULL));

  ASSERT_CALL(sleqp_scaling_create(&scaling,
                                   constrained_num_variables,
                                   constrained_num_constraints));

  ASSERT_CALL(
    sleqp_func_scaling_from_gradient(scaling,
                                     sleqp_iterate_get_func_grad(iterate),
                                     eps));

  ASSERT_CALL(sleqp_scaling_from_cons_jac(scaling,
                                          sleqp_iterate_get_cons_jac(iterate),
                                          eps));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  params,
                                  options,
                                  constrained_initial,
                                  scaling));

  // 100 iterations should be plenty...
  ASSERT_CALL(sleqp_solver_solve(solver, 100, -1));

  SleqpIterate* solution_iterate;

  ASSERT_CALL(sleqp_solver_get_solution(solver, &solution_iterate));

  ck_assert_int_eq(sleqp_solver_get_status(solver), SLEQP_STATUS_OPTIMAL);

  SleqpSparseVec* actual_solution = sleqp_iterate_get_primal(solution_iterate);

  ck_assert(sleqp_sparse_vector_eq(actual_solution, constrained_optimal, 1e-6));

  ASSERT_CALL(sleqp_solver_release(&solver));

  ASSERT_CALL(sleqp_scaling_release(&scaling));

  ASSERT_CALL(sleqp_iterate_release(&iterate));
}
END_TEST

Suite*
constrained_test_suite()
{
  Suite* suite;
  TCase* tc_cons;

  suite = suite_create("Constrained tests");

  tc_cons = tcase_create("Constrained solution test");

  tcase_add_checked_fixture(tc_cons,
                            constrained_test_setup,
                            constrained_test_teardown);

  tcase_add_test(tc_cons, test_solve);

  tcase_add_test(tc_cons, test_exact_linesearch);

  tcase_add_test(tc_cons, test_parametric_solve);

  tcase_add_test(tc_cons, test_sr1_solve);

  tcase_add_test(tc_cons, test_bfgs_solve_no_sizing);

  tcase_add_test(tc_cons, test_bfgs_solve_centered_ol_sizing);

  tcase_add_test(tc_cons, test_unscaled_solve);

  tcase_add_test(tc_cons, test_scaled_solve);

  tcase_add_test(tc_cons, test_scaled_sr1_solve);

  tcase_add_test(tc_cons, test_scaled_bfgs_solve);

  tcase_add_test(tc_cons, test_auto_scaled_solve);

  suite_add_tcase(suite, tc_cons);

  return suite;
}

TEST_MAIN(constrained_test_suite)
