#include <check.h>
#include <stdlib.h>

#include "cmp.h"
#include "mem.h"
#include "solver.h"

#include "test_common.h"

#include "rosenbrock_fixture.h"

SleqpSettings* settings;
SleqpProblem* problem;
SleqpSolver* solver;

void
callback_setup()
{
  rosenbrock_setup();

  ASSERT_CALL(sleqp_settings_create(&settings));

  ASSERT_CALL(sleqp_problem_create_simple(&problem,
                                          rosenbrock_func,
                                          rosenbrock_var_lb,
                                          rosenbrock_var_ub,
                                          rosenbrock_cons_lb,
                                          rosenbrock_cons_ub,
                                          settings));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  rosenbrock_initial,
                                  NULL));
}

void
callback_teardown()
{
  ASSERT_CALL(sleqp_solver_release(&solver));

  ASSERT_CALL(sleqp_problem_release(&problem));

  ASSERT_CALL(sleqp_settings_release(&settings));

  rosenbrock_teardown();
}

typedef struct CallbackData
{
  bool called;
  bool accepted_optimal_sol;
} CallbackData;

SLEQP_RETCODE
accepted_iterate(SleqpSolver* solver, SleqpIterate* iterate, void* data)
{
  CallbackData* callback_data = (CallbackData*)data;

  callback_data->called = true;

  SleqpVec* solution = sleqp_iterate_primal(iterate);

  assert(solution->dim == 2);

  if (sleqp_vec_eq(solution, rosenbrock_optimum, 1e-6))
  {
    callback_data->accepted_optimal_sol = true;
  }

  return SLEQP_OKAY;
}

START_TEST(test_add_invalid_number)
{
  sleqp_log_set_level(SLEQP_LOG_SILENT);

  ck_assert_int_eq(sleqp_solver_add_callback(solver, -1, NULL, NULL),
                   SLEQP_ERROR);

  ck_assert_int_eq(sleqp_error_type(), SLEQP_ILLEGAL_ARGUMENT);
}
END_TEST

START_TEST(test_remove_invalid_number)
{
  sleqp_log_set_level(SLEQP_LOG_SILENT);

  ck_assert_int_eq(sleqp_solver_remove_callback(solver, -1, NULL, NULL),
                   SLEQP_ERROR);

  ck_assert_int_eq(sleqp_error_type(), SLEQP_ILLEGAL_ARGUMENT);
}
END_TEST

START_TEST(test_simple_callback)
{
  CallbackData callback_data;

  callback_data.called               = false;
  callback_data.accepted_optimal_sol = false;

  ASSERT_CALL(sleqp_solver_add_callback(solver,
                                        SLEQP_SOLVER_EVENT_ACCEPTED_ITERATE,
                                        (void*)accepted_iterate,
                                        &callback_data));

  // 100 iterations should be plenty...
  ASSERT_CALL(sleqp_solver_solve(solver, 100, -1));

  ck_assert_int_eq(callback_data.called, true);
  ck_assert_int_eq(callback_data.accepted_optimal_sol, true);
}
END_TEST

START_TEST(test_remove_callback)
{
  CallbackData callback_data;

  callback_data.called               = false;
  callback_data.accepted_optimal_sol = false;

  ASSERT_CALL(sleqp_solver_add_callback(solver,
                                        SLEQP_SOLVER_EVENT_ACCEPTED_ITERATE,
                                        accepted_iterate,
                                        &callback_data));

  ASSERT_CALL(sleqp_solver_remove_callback(solver,
                                           SLEQP_SOLVER_EVENT_ACCEPTED_ITERATE,
                                           accepted_iterate,
                                           &callback_data));

  // 100 iterations should be plenty...
  ASSERT_CALL(sleqp_solver_solve(solver, 100, -1));

  ck_assert_int_eq(callback_data.called, false);
  ck_assert_int_eq(callback_data.accepted_optimal_sol, false);
}
END_TEST

Suite*
callback_test_suite()
{
  Suite* suite;
  TCase* tc_callback;

  suite = suite_create("Callback tests");

  tc_callback = tcase_create("Callback test");

  tcase_add_checked_fixture(tc_callback, callback_setup, callback_teardown);

  tcase_add_test(tc_callback, test_simple_callback);
  tcase_add_test(tc_callback, test_remove_callback);

  tcase_add_test(tc_callback, test_add_invalid_number);
  tcase_add_test(tc_callback, test_remove_invalid_number);

  suite_add_tcase(suite, tc_callback);

  return suite;
}

TEST_MAIN(callback_test_suite)
