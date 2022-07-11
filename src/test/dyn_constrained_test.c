#include <check.h>
#include <stdlib.h>

#include "cmp.h"
#include "constrained_fixture.h"
#include "mem.h"
#include "scale.h"
#include "solver.h"

#include "test_common.h"

#include "dyn_constrained_fixture.h"

SleqpParams* params;
SleqpOptions* options;
SleqpProblem* problem;

void
setup()
{
  dyn_constrained_setup();

  ASSERT_CALL(sleqp_params_create(&params));

  ASSERT_CALL(sleqp_options_create(&options));

  ASSERT_CALL(sleqp_problem_create_simple(&problem,
                                          dyn_constrained_func,
                                          params,
                                          constrained_var_lb,
                                          constrained_var_ub,
                                          constrained_cons_lb,
                                          constrained_cons_ub));
}

void
teardown()
{
  ASSERT_CALL(sleqp_problem_release(&problem));

  ASSERT_CALL(sleqp_options_release(&options));

  ASSERT_CALL(sleqp_params_release(&params));

  dyn_constrained_teardown();
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

  ck_assert(sleqp_vec_eq(actual_solution, constrained_optimum, 1e-6));

  ASSERT_CALL(sleqp_solver_release(&solver));
}

START_TEST(test_solve)
{
  SleqpSolver* solver;

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  params,
                                  options,
                                  constrained_initial,
                                  NULL));

  solve_and_release_solver(solver);
}
END_TEST

START_TEST(test_scaled_solve)
{
  SleqpSolver* solver;

  SleqpScaling* scaling;

  ASSERT_CALL(sleqp_scaling_create(&scaling,
                                   constrained_num_variables,
                                   constrained_num_constraints));

  ASSERT_CALL(sleqp_scaling_set_obj_weight(scaling, 2));

  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling, 0, -5));
  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling, 1, 5));

  ASSERT_CALL(sleqp_scaling_set_cons_weight(scaling, 0, -1));
  ASSERT_CALL(sleqp_scaling_set_cons_weight(scaling, 1, -2));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  params,
                                  options,
                                  constrained_initial,
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

  suite = suite_create("Dynamic constrained tests");

  tc_dyn = tcase_create("Dynamic constrained solution test");

  tcase_add_checked_fixture(tc_dyn, setup, teardown);

  tcase_add_test(tc_dyn, test_solve);
  tcase_add_test(tc_dyn, test_scaled_solve);

  suite_add_tcase(suite, tc_dyn);

  return suite;
}

TEST_MAIN(dyn_test_suite)
