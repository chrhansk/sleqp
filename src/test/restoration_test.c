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

SleqpParams* params;

SleqpProblem* problem;
SleqpIterate* iterate;

SleqpProblem* restoration_problem;

void
restoration_setup()
{
  wachbieg_setup();

  ASSERT_CALL(sleqp_params_create(&params));

  ASSERT_CALL(sleqp_problem_create_simple(&problem,
                                          wachbieg_func,
                                          params,
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
    sleqp_restoration_problem_create(&restoration_problem, params, problem));
}

void
restoration_teardown()
{
  ASSERT_CALL(sleqp_problem_release(&restoration_problem));

  ASSERT_CALL(sleqp_iterate_release(&iterate));

  ASSERT_CALL(sleqp_problem_release(&problem));

  ASSERT_CALL(sleqp_params_release(&params));

  wachbieg_teardown();
}

START_TEST(test_transform)
{
  SleqpSparseVec* restoration_primal;

  SleqpSparseVec* residuals;

  ASSERT_CALL(sleqp_sparse_vector_create_empty(
    &restoration_primal,
    sleqp_problem_num_vars(restoration_problem)));

  ASSERT_CALL(
    sleqp_sparse_vector_create_empty(&residuals,
                                     sleqp_problem_num_cons(problem)));

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

  const double eps = sleqp_params_value(params, SLEQP_PARAM_EPS);

  ck_assert(sleqp_is_eq(sleqp_iterate_obj_val(restoration_iterate),
                        .5 * sleqp_sparse_vector_norm_sq(residuals),
                        eps));

  ASSERT_CALL(sleqp_iterate_release(&restoration_iterate));

  ASSERT_CALL(sleqp_sparse_vector_free(&residuals));

  ASSERT_CALL(sleqp_sparse_vector_free(&restoration_primal));
}
END_TEST

START_TEST(test_restore)
{
  const double zero_eps = sleqp_params_value(params, SLEQP_PARAM_ZERO_EPS);

  SleqpSparseVec* transformed_primal;
  SleqpSparseVec* restored_primal;

  ASSERT_CALL(sleqp_sparse_vector_create_empty(
    &transformed_primal,
    sleqp_problem_num_vars(restoration_problem)));

  ASSERT_CALL(
    sleqp_sparse_vector_create_empty(&restored_primal,
                                     sleqp_problem_num_vars(problem)));

  ASSERT_CALL(
    sleqp_restoration_problem_transform(problem,
                                        sleqp_iterate_primal(iterate),
                                        sleqp_iterate_cons_val(iterate),
                                        transformed_primal));

  ASSERT_CALL(sleqp_restoration_problem_restore(problem,
                                                transformed_primal,
                                                restored_primal));

  ck_assert(sleqp_sparse_vector_eq(sleqp_iterate_primal(iterate),
                                   restored_primal,
                                   zero_eps));

  ASSERT_CALL(sleqp_sparse_vector_free(&restored_primal));

  ASSERT_CALL(sleqp_sparse_vector_free(&transformed_primal));
}
END_TEST

START_TEST(test_solve)
{
  SleqpSparseVec* initial;
  SleqpSparseVec* residuals;
  SleqpOptions* options;
  SleqpSolver* solver;

  ASSERT_CALL(sleqp_sparse_vector_create_empty(
    &initial,
    sleqp_problem_num_vars(restoration_problem)));

  ASSERT_CALL(
    sleqp_sparse_vector_create_empty(&residuals,
                                     sleqp_problem_num_cons(problem)));

  ASSERT_CALL(sleqp_options_create(&options));

  ASSERT_CALL(sleqp_options_set_enum_value(options,
                                           SLEQP_OPTION_ENUM_DERIV_CHECK,
                                           SLEQP_DERIV_CHECK_FIRST));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  restoration_problem,
                                  params,
                                  options,
                                  initial,
                                  NULL));

  // 100 iterations should be plenty...
  ASSERT_CALL(sleqp_solver_solve(solver, 100, -1));

  SleqpIterate* solution_iterate;

  ASSERT_CALL(sleqp_solver_solution(solver, &solution_iterate));

  ck_assert_int_eq(sleqp_solver_status(solver), SLEQP_STATUS_OPTIMAL);

  SleqpSparseVec* transformed_solution = sleqp_iterate_primal(solution_iterate);

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

  ck_assert(sleqp_is_zero(sleqp_sparse_vector_inf_norm(residuals), 1e-4));

  ASSERT_CALL(sleqp_solver_release(&solver));

  ASSERT_CALL(sleqp_options_release(&options));

  ASSERT_CALL(sleqp_sparse_vector_free(&residuals));

  ASSERT_CALL(sleqp_sparse_vector_free(&initial));
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
