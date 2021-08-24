#include <stdlib.h>
#include <check.h>

#include "cmp.h"
#include "mem.h"
#include "solver.h"

#include "test_common.h"

#include "log_rosenbrock_fixture.h"

void test_step_rule(SLEQP_STEP_RULE step_rule)
{
  SleqpParams* params;
  SleqpOptions* options;
  SleqpProblem* problem;
  SleqpSolver* solver;

  ASSERT_CALL(sleqp_params_create(&params));

  ASSERT_CALL(sleqp_options_create(&options));

  ASSERT_CALL(sleqp_options_set_int(options,
                                    SLEQP_OPTION_INT_STEP_RULE,
                                    step_rule));

  ASSERT_CALL(sleqp_problem_create_simple(&problem,
                                          log_rosenbrock_func,
                                          params,
                                          log_rosenbrock_var_lb,
                                          log_rosenbrock_var_ub,
                                          log_rosenbrock_cons_lb,
                                          log_rosenbrock_cons_ub));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  params,
                                  options,
                                  log_rosenbrock_initial,
                                  NULL));

  ASSERT_CALL(sleqp_solver_solve(solver, 1000, -1));

  SleqpIterate* solution_iterate;

  ASSERT_CALL(sleqp_solver_get_solution(solver,
                                        &solution_iterate));

  ck_assert_int_eq(sleqp_solver_get_status(solver), SLEQP_STATUS_OPTIMAL);

  SleqpSparseVec* actual_solution = sleqp_iterate_get_primal(solution_iterate);

  ck_assert(sleqp_sparse_vector_eq(actual_solution,
                                   log_rosenbrock_optimal,
                                   1e-6));

  ASSERT_CALL(sleqp_solver_release(&solver));

  ASSERT_CALL(sleqp_problem_release(&problem));

  ASSERT_CALL(sleqp_options_release(&options));

  ASSERT_CALL(sleqp_params_release(&params));
}

START_TEST(test_direct)
{
  test_step_rule(SLEQP_STEP_RULE_DIRECT);
}
END_TEST

START_TEST(test_window)
{
  test_step_rule(SLEQP_STEP_RULE_WINDOW);
}
END_TEST

START_TEST(test_minstep)
{
  test_step_rule(SLEQP_STEP_RULE_MINSTEP);
}
END_TEST

Suite* unconstrained_test_suite()
{
  Suite *suite;
  TCase *tc_step_rule;

  suite = suite_create("Unconstrained tests");

  tc_step_rule = tcase_create("Unconstrained solution test");

  tcase_add_checked_fixture(tc_step_rule,
                            log_rosenbrock_setup,
                            log_rosenbrock_teardown);

  tcase_add_test(tc_step_rule, test_direct);
  tcase_add_test(tc_step_rule, test_window);
  tcase_add_test(tc_step_rule, test_minstep);

  suite_add_tcase(suite, tc_step_rule);

  return suite;
}

TEST_MAIN(unconstrained_test_suite)
