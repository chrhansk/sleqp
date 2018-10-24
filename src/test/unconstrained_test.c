#include <stdlib.h>
#include <check.h>

#include "sleqp_cmp.h"
#include "sleqp_mem.h"
#include "sleqp_solver.h"

#include "test_common.h"

#include "rosenbrock_fixture.h"


START_TEST(test_unconstrained_cauchy_direction)
{
  SleqpProblem* problem;
  SleqpSolver* solver;

  ASSERT_CALL(sleqp_problem_create(&problem,
                                   rosenbrock_func,
                                   rosenbrock_var_lb,
                                   rosenbrock_var_ub,
                                   rosenbrock_cons_lb,
                                   rosenbrock_cons_ub));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  rosenbrock_x));

  ASSERT_CALL(sleqp_solve(solver));

  ASSERT_CALL(sleqp_solver_free(&solver));

  ASSERT_CALL(sleqp_problem_free(&problem));
}
END_TEST

Suite* unconstrained_test_suite()
{
  Suite *suite;
  TCase *tc_uncons;

  suite = suite_create("Unconstrained tests");

  tc_uncons = tcase_create("Unconstrained Cauchy direction");

  tcase_add_checked_fixture(tc_uncons,
                            rosenbrock_setup,
                            rosenbrock_teardown);

  tcase_add_test(tc_uncons, test_unconstrained_cauchy_direction);
  suite_add_tcase(suite, tc_uncons);

  return suite;
}


int main()
{
  int num_fails;
  Suite* suite;
  SRunner* srunner;

  suite = unconstrained_test_suite();
  srunner = srunner_create(suite);

  srunner_set_fork_status(srunner, CK_NOFORK);
  srunner_run_all(srunner, CK_NORMAL);

  num_fails = srunner_ntests_failed(srunner);

  srunner_free(srunner);

  return (num_fails > 0) ? EXIT_FAILURE : EXIT_SUCCESS;
}
