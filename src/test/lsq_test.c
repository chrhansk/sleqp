#include <stdlib.h>
#include <check.h>
#include <math.h>

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

void lsq_setup()
{
  rosenbrock_setup();

  rosenbrock_lsq_setup();

  ASSERT_CALL(sleqp_params_create(&params));
  ASSERT_CALL(sleqp_options_create(&options));
}

void lsq_teardown()
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

  ASSERT_CALL(sleqp_options_set_int(options,
                                    SLEQP_OPTION_INT_DERIV_CHECK,
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

  ASSERT_CALL(sleqp_solver_get_solution(solver,
                                        &solution_iterate));

  ck_assert_int_eq(sleqp_solver_get_status(solver), SLEQP_STATUS_OPTIMAL);

  SleqpSparseVec* actual_solution = sleqp_iterate_get_primal(solution_iterate);

  ck_assert(sleqp_sparse_vector_eq(actual_solution,
                                   rosenbrock_optimal,
                                   1e-6));

  ASSERT_CALL(sleqp_solver_release(&solver));

  ASSERT_CALL(sleqp_problem_release(&problem));
}
END_TEST

START_TEST(test_lsqr_solve)
{
  SleqpProblem* problem;
  SleqpSolver* solver;

  ASSERT_CALL(sleqp_options_set_int(options,
                                    SLEQP_OPTION_INT_TR_SOLVER,
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

  ASSERT_CALL(sleqp_solver_get_solution(solver,
                                        &solution_iterate));

  ck_assert_int_eq(sleqp_solver_get_status(solver), SLEQP_STATUS_OPTIMAL);

  SleqpSparseVec* actual_solution = sleqp_iterate_get_primal(solution_iterate);

  ck_assert(sleqp_sparse_vector_eq(actual_solution,
                                   rosenbrock_optimal,
                                   1e-6));

  ASSERT_CALL(sleqp_solver_release(&solver));

  ASSERT_CALL(sleqp_problem_release(&problem));
}
END_TEST

START_TEST(test_constrained_lsqr_solve)
{
  SleqpProblem* problem;
  SleqpSolver* solver;

  SleqpSparseMatrix* linear_coeffs;
  SleqpSparseVec* linear_lb;
  SleqpSparseVec* linear_ub;

  const double inf = sleqp_infinity();

  ASSERT_CALL(sleqp_sparse_matrix_create(&linear_coeffs, 1, rosenbrock_num_variables, 2));

  ASSERT_CALL(sleqp_sparse_matrix_push_column(linear_coeffs, 0));

  ASSERT_CALL(sleqp_sparse_matrix_push(linear_coeffs, 0, 0, 1.));

  ASSERT_CALL(sleqp_sparse_matrix_push_column(linear_coeffs, 1));

  ASSERT_CALL(sleqp_sparse_matrix_push(linear_coeffs, 0, 1, 1.));

  ASSERT_CALL(sleqp_sparse_vector_create(&linear_lb, 1, 1));
  ASSERT_CALL(sleqp_sparse_vector_push(linear_lb, 0, 1.));

  ASSERT_CALL(sleqp_sparse_vector_create(&linear_ub, 1, 1));
  ASSERT_CALL(sleqp_sparse_vector_push(linear_ub, 0, inf));

  ASSERT_CALL(sleqp_options_set_int(options,
                                    SLEQP_OPTION_INT_TR_SOLVER,
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

  ASSERT_CALL(sleqp_sparse_vector_clear(rosenbrock_initial));
  ASSERT_CALL(sleqp_sparse_vector_reserve(rosenbrock_initial, 2));

  ASSERT_CALL(sleqp_sparse_vector_push(rosenbrock_initial, 0, -1.));
  ASSERT_CALL(sleqp_sparse_vector_push(rosenbrock_initial, 1, -1.));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  params,
                                  options,
                                  rosenbrock_initial,
                                  NULL));

  // 100 iterations should be plenty...
  ASSERT_CALL(sleqp_solver_solve(solver, 100, -1));

  SleqpIterate* solution_iterate;

  ASSERT_CALL(sleqp_solver_get_solution(solver,
                                        &solution_iterate));

  ck_assert_int_eq(sleqp_solver_get_status(solver), SLEQP_STATUS_OPTIMAL);

  SleqpSparseVec* actual_solution = sleqp_iterate_get_primal(solution_iterate);

  ck_assert(sleqp_sparse_vector_eq(actual_solution,
                                   rosenbrock_optimal,
                                   1e-6));

  ASSERT_CALL(sleqp_solver_release(&solver));

  ASSERT_CALL(sleqp_problem_release(&problem));

  ASSERT_CALL(sleqp_sparse_vector_free(&linear_ub));
  ASSERT_CALL(sleqp_sparse_vector_free(&linear_lb));

  ASSERT_CALL(sleqp_sparse_matrix_release(&linear_coeffs));
}
END_TEST

START_TEST(test_scaled_solve)
{
  SleqpProblem* problem;
  SleqpScaling* scaling;
  SleqpSolver* solver;

  ASSERT_CALL(sleqp_scaling_create(&scaling,
                                   rosenbrock_num_variables,
                                   rosenbrock_num_constraints));

  ASSERT_CALL(sleqp_scaling_set_func_weight(scaling,
                                            -2));

  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling,
                                           0,
                                           -1));

  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling,
                                           1,
                                           -1));

  ASSERT_CALL(sleqp_options_set_int(options,
                                    SLEQP_OPTION_INT_DERIV_CHECK,
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

  ASSERT_CALL(sleqp_solver_get_solution(solver,
                                        &solution_iterate));

  ck_assert_int_eq(sleqp_solver_get_status(solver), SLEQP_STATUS_OPTIMAL);

  SleqpSparseVec* actual_solution = sleqp_iterate_get_primal(solution_iterate);

  ck_assert(sleqp_sparse_vector_eq(actual_solution,
                                   rosenbrock_optimal,
                                   1e-6));

  ASSERT_CALL(sleqp_solver_release(&solver));

  ASSERT_CALL(sleqp_problem_release(&problem));

  ASSERT_CALL(sleqp_scaling_release(&scaling));
}
END_TEST


Suite* lsq_test_suite()
{
  Suite *suite;
  TCase *tc_uncons;

  suite = suite_create("LSQ tests");

  tc_uncons = tcase_create("LSQ solution test");

  tcase_add_checked_fixture(tc_uncons,
                            lsq_setup,
                            lsq_teardown);

  tcase_add_test(tc_uncons, test_unconstrained_solve);
  tcase_add_test(tc_uncons, test_lsqr_solve);
  tcase_add_test(tc_uncons, test_constrained_lsqr_solve);
  tcase_add_test(tc_uncons, test_scaled_solve);

  suite_add_tcase(suite, tc_uncons);

  return suite;
}

TEST_MAIN(lsq_test_suite)
