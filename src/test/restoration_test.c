#include <check.h>
#include <stdlib.h>

#include "cmp.h"
#include "feas.h"
#include "mem.h"
#include "restoration.h"
#include "solver.h"
#include "util.h"

#include "test_common.h"

#include "wachbieg_fixture.h"

SleqpSettings* settings;

SleqpProblem* problem;
SleqpIterate* iterate;

SleqpProblem* restoration_problem;

void
restoration_setup()
{
  wachbieg_setup();

  ASSERT_CALL(sleqp_settings_create(&settings));

  ASSERT_CALL(sleqp_problem_create_simple(&problem,
                                          wachbieg_func,
                                          settings,
                                          wachbieg_var_lb,
                                          wachbieg_var_ub,
                                          wachbieg_cons_lb,
                                          wachbieg_cons_ub));

  ASSERT_CALL(sleqp_iterate_create(&iterate, problem, wachbieg_initial));

  bool reject;

  ASSERT_CALL(
    sleqp_set_and_evaluate(problem, iterate, SLEQP_VALUE_REASON_NONE, &reject));

  assert(!reject);

  ASSERT_CALL(
    sleqp_restoration_problem_create(&restoration_problem, settings, problem));
}

void
restoration_teardown()
{
  ASSERT_CALL(sleqp_problem_release(&restoration_problem));

  ASSERT_CALL(sleqp_iterate_release(&iterate));

  ASSERT_CALL(sleqp_problem_release(&problem));

  ASSERT_CALL(sleqp_settings_release(&settings));

  wachbieg_teardown();
}

START_TEST(test_transform)
{
  SleqpVec* restoration_primal;

  SleqpVec* residuals;

  ASSERT_CALL(
    sleqp_vec_create_empty(&restoration_primal,
                           sleqp_problem_num_vars(restoration_problem)));

  ASSERT_CALL(
    sleqp_vec_create_empty(&residuals, sleqp_problem_num_cons(problem)));

  ASSERT_CALL(sleqp_feasibility_residuals(problem,
                                          sleqp_iterate_cons_val(iterate),
                                          residuals,
                                          NULL));

  ASSERT_CALL(
    sleqp_restoration_problem_transform(problem,
                                        sleqp_iterate_primal(iterate),
                                        sleqp_iterate_cons_val(iterate),
                                        restoration_primal));

  SleqpIterate* restoration_iterate;

  ASSERT_CALL(sleqp_iterate_create(&restoration_iterate,
                                   restoration_problem,
                                   restoration_primal));

  bool reject;

  ASSERT_CALL(sleqp_set_and_evaluate(restoration_problem,
                                     restoration_iterate,
                                     SLEQP_VALUE_REASON_NONE,
                                     &reject));

  assert(!reject);

  const double eps = sleqp_settings_real_value(settings, SLEQP_SETTINGS_REAL_EPS);

  ck_assert(sleqp_is_eq(sleqp_iterate_obj_val(restoration_iterate),
                        .5 * sleqp_vec_norm_sq(residuals),
                        eps));

  ASSERT_CALL(sleqp_iterate_release(&restoration_iterate));

  ASSERT_CALL(sleqp_vec_free(&residuals));

  ASSERT_CALL(sleqp_vec_free(&restoration_primal));
}
END_TEST

START_TEST(test_restore)
{
  const double zero_eps = sleqp_settings_real_value(settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  SleqpVec* transformed_primal;
  SleqpVec* restored_primal;

  ASSERT_CALL(
    sleqp_vec_create_empty(&transformed_primal,
                           sleqp_problem_num_vars(restoration_problem)));

  ASSERT_CALL(
    sleqp_vec_create_empty(&restored_primal, sleqp_problem_num_vars(problem)));

  ASSERT_CALL(
    sleqp_restoration_problem_transform(problem,
                                        sleqp_iterate_primal(iterate),
                                        sleqp_iterate_cons_val(iterate),
                                        transformed_primal));

  ASSERT_CALL(sleqp_restoration_problem_restore(problem,
                                                transformed_primal,
                                                restored_primal));

  ck_assert(
    sleqp_vec_eq(sleqp_iterate_primal(iterate), restored_primal, zero_eps));

  ASSERT_CALL(sleqp_vec_free(&restored_primal));

  ASSERT_CALL(sleqp_vec_free(&transformed_primal));
}
END_TEST

START_TEST(test_solve)
{
  SleqpVec* initial;
  SleqpVec* residuals;
  SleqpSettings* settings;
  SleqpSolver* solver;

  ASSERT_CALL(
    sleqp_vec_create_empty(&initial,
                           sleqp_problem_num_vars(restoration_problem)));

  ASSERT_CALL(
    sleqp_vec_create_empty(&residuals, sleqp_problem_num_cons(problem)));

  ASSERT_CALL(sleqp_settings_create(&settings));

  ASSERT_CALL(sleqp_settings_set_enum_value(settings,
                                           SLEQP_SETTINGS_ENUM_DERIV_CHECK,
                                           SLEQP_DERIV_CHECK_FIRST));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  restoration_problem,
                                  settings,
                                  initial,
                                  NULL));

  // 100 iterations should be plenty...
  ASSERT_CALL(sleqp_solver_solve(solver, 100, -1));

  SleqpIterate* solution_iterate;

  ASSERT_CALL(sleqp_solver_solution(solver, &solution_iterate));

  ck_assert_int_eq(sleqp_solver_status(solver), SLEQP_STATUS_OPTIMAL);

  SleqpVec* transformed_solution = sleqp_iterate_primal(solution_iterate);

  ASSERT_CALL(sleqp_restoration_problem_restore(problem,
                                                transformed_solution,
                                                sleqp_iterate_primal(iterate)));

  bool reject;

  ASSERT_CALL(
    sleqp_set_and_evaluate(problem, iterate, SLEQP_VALUE_REASON_NONE, &reject));

  assert(!reject);

  ASSERT_CALL(sleqp_feasibility_residuals(problem,
                                          sleqp_iterate_cons_val(iterate),
                                          residuals,
                                          NULL));

  ck_assert(sleqp_is_zero(sleqp_vec_inf_norm(residuals), 1e-4));

  ASSERT_CALL(sleqp_solver_release(&solver));

  ASSERT_CALL(sleqp_settings_release(&settings));

  ASSERT_CALL(sleqp_vec_free(&residuals));

  ASSERT_CALL(sleqp_vec_free(&initial));
}
END_TEST

Suite*
restoration_test_suite()
{
  Suite* suite;
  TCase* tc_restoration;

  suite = suite_create("Restoration tests");

  tc_restoration = tcase_create("Restoration test");

  tcase_add_checked_fixture(tc_restoration,
                            restoration_setup,
                            restoration_teardown);

  tcase_add_test(tc_restoration, test_transform);

  tcase_add_test(tc_restoration, test_restore);

  tcase_add_test(tc_restoration, test_solve);

  suite_add_tcase(suite, tc_restoration);

  return suite;
}

TEST_MAIN(restoration_test_suite)
