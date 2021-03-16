#include <stdlib.h>
#include <check.h>

#include "sleqp_cmp.h"
#include "sleqp_mem.h"
#include "sleqp_solver.h"

#include "test_common.h"

#include "rosenbrock_fixture.h"

SleqpParams* params;
SleqpOptions* options;
SleqpProblem* problem;
SleqpSolver* solver;

void callback_setup()
{
  rosenbrock_setup();

  ASSERT_CALL(sleqp_params_create(&params));

  ASSERT_CALL(sleqp_options_create(&options));

  ASSERT_CALL(sleqp_problem_create(&problem,
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
}

void callback_teardown()
{
  ASSERT_CALL(sleqp_solver_release(&solver));

  ASSERT_CALL(sleqp_problem_free(&problem));

  ASSERT_CALL(sleqp_options_release(&options));

  ASSERT_CALL(sleqp_params_release(&params));

  rosenbrock_teardown();
}

typedef struct CallbackData
{
  bool called;
  bool accepted_optimal_sol;
} CallbackData;

SLEQP_RETCODE accepted_iterate(SleqpSolver* solver,
                               SleqpIterate* iterate,
                               void* data)
{
  CallbackData* callback_data = (CallbackData*) data;

  callback_data->called = true;

  SleqpSparseVec* solution = sleqp_iterate_get_primal(iterate);

  assert(solution->dim == 2);

  if(sleqp_sparse_vector_eq(solution, rosenbrock_optimal, 1e-6))
  {
    callback_data->accepted_optimal_sol = true;
  }

  return SLEQP_OKAY;
}

START_TEST(test_simple_callback)
{
  CallbackData callback_data;

  callback_data.called = false;
  callback_data.accepted_optimal_sol = false;

  ASSERT_CALL(sleqp_solver_add_callback(solver,
                                        SLEQP_SOLVER_EVENT_ACCEPTED_ITERATE,
                                        (void*) accepted_iterate,
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

  callback_data.called = false;
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

Suite* callback_test_suite()
{
  Suite *suite;
  TCase *tc_callback;

  suite = suite_create("Callback tests");

  tc_callback = tcase_create("Callback test");

  tcase_add_checked_fixture(tc_callback,
                            callback_setup,
                            callback_teardown);

  tcase_add_test(tc_callback, test_simple_callback);
  tcase_add_test(tc_callback, test_remove_callback);

  suite_add_tcase(suite, tc_callback);

  return suite;
}

TEST_MAIN(callback_test_suite)
