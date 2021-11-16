#include <check.h>
#include <stdlib.h>

#include "cmp.h"
#include "mem.h"
#include "solver.h"

#include "test_common.h"

#include "rosenbrock_fixture.h"

START_TEST(test_unconstrained_solve)
{
  SleqpParams* params;
  SleqpOptions* options;
  SleqpProblem* problem;
  SleqpSolver* solver;

  ASSERT_CALL(sleqp_params_create(&params));

  ASSERT_CALL(sleqp_options_create(&options));

  ASSERT_CALL(sleqp_problem_create_simple(&problem,
                                          rosenbrock_func,
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

  SleqpSparseVec* actual_solution = sleqp_iterate_primal(solution_iterate);

  ck_assert(sleqp_sparse_vector_eq(actual_solution, rosenbrock_optimal, 1e-6));

  ASSERT_CALL(sleqp_solver_release(&solver));

  ASSERT_CALL(sleqp_problem_release(&problem));

  ASSERT_CALL(sleqp_options_release(&options));

  ASSERT_CALL(sleqp_params_release(&params));
}
END_TEST

Suite*
unconstrained_test_suite()
{
  Suite* suite;
  TCase* tc_uncons;

  suite = suite_create("Unconstrained tests");

  tc_uncons = tcase_create("Unconstrained solution test");

  tcase_add_checked_fixture(tc_uncons, rosenbrock_setup, rosenbrock_teardown);

  tcase_add_test(tc_uncons, test_unconstrained_solve);
  suite_add_tcase(suite, tc_uncons);

  return suite;
}

TEST_MAIN(unconstrained_test_suite)
