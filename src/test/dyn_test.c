#include <check.h>
#include <stdlib.h>

#include "cmp.h"
#include "mem.h"
#include "scale.h"
#include "solver.h"

#include "test_common.h"

#include "dyn_rosenbrock_fixture.h"

SleqpSettings* settings;
SleqpProblem* problem;

void
setup()
{
  dyn_rosenbrock_setup();

  ASSERT_CALL(sleqp_settings_create(&settings));

  ASSERT_CALL(sleqp_problem_create_simple(&problem,
                                          dyn_rosenbrock_func,
                                          settings,
                                          rosenbrock_var_lb,
                                          rosenbrock_var_ub,
                                          rosenbrock_cons_lb,
                                          rosenbrock_cons_ub));
}

void
teardown()
{
  ASSERT_CALL(sleqp_problem_release(&problem));

  ASSERT_CALL(sleqp_settings_release(&settings));

  dyn_rosenbrock_teardown();
}

void
solve_and_release_solver(SleqpSolver* solver)
{
  // 1000 iterations, one minute time limit
  ASSERT_CALL(sleqp_solver_solve(solver, 1000, 60.));

  SleqpIterate* iterate;

  ASSERT_CALL(sleqp_solver_solution(solver, &iterate));

  ck_assert_int_eq(sleqp_solver_status(solver), SLEQP_STATUS_OPTIMAL);

  SleqpVec* actual_solution = sleqp_iterate_primal(iterate);

  ck_assert(sleqp_vec_eq(actual_solution, rosenbrock_optimum, 1e-6));

  ASSERT_CALL(sleqp_solver_release(&solver));
}

START_TEST(test_solve)
{
  SleqpSolver* solver;

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  settings,
                                  rosenbrock_initial,
                                  NULL));

  solve_and_release_solver(solver);
}
END_TEST

START_TEST(test_scaled_solve)
{
  SleqpSolver* solver;

  SleqpScaling* scaling;

  ASSERT_CALL(
    sleqp_scaling_create(&scaling, rosenbrock_num_vars, rosenbrock_num_cons));

  ASSERT_CALL(sleqp_scaling_set_obj_weight(scaling, 2));

  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling, 0, -5));
  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling, 1, 5));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  settings,
                                  rosenbrock_initial,
                                  scaling));

  solve_and_release_solver(solver);

  ASSERT_CALL(sleqp_scaling_release(&scaling));
}
END_TEST

Suite*
dyn_test_suite()
{
  Suite* suite;
  TCase* tc_dyn;

  suite = suite_create("Dynamic unconstrained tests");

  tc_dyn = tcase_create("Dynamic unconstrained solution test");

  tcase_add_checked_fixture(tc_dyn, setup, teardown);

  tcase_add_test(tc_dyn, test_solve);
  tcase_add_test(tc_dyn, test_scaled_solve);

  suite_add_tcase(suite, tc_dyn);

  return suite;
}

TEST_MAIN(dyn_test_suite)
