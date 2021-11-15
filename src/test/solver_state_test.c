#include <check.h>

#include "constrained_fixture.h"
#include "test_common.h"

#include "solver.h"

SleqpParams* params;
SleqpOptions* options;
SleqpProblem* problem;
SleqpSolver* solver;

void
solver_state_setup()
{
  constrained_setup();

  ASSERT_CALL(sleqp_params_create(&params));

  ASSERT_CALL(sleqp_options_create(&options));

  ASSERT_CALL(sleqp_problem_create_simple(&problem,
                                          constrained_func,
                                          params,
                                          constrained_var_lb,
                                          constrained_var_ub,
                                          constrained_cons_lb,
                                          constrained_cons_ub));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  params,
                                  options,
                                  constrained_initial,
                                  NULL));
}

void
solver_state_teardown()
{
  ASSERT_CALL(sleqp_solver_release(&solver));

  ASSERT_CALL(sleqp_problem_release(&problem));

  ASSERT_CALL(sleqp_options_release(&options));

  ASSERT_CALL(sleqp_params_release(&params));

  constrained_teardown();
}

START_TEST(test_stationarity_residuals)
{
  SleqpSparseVec* stationarity_residuals;

  ASSERT_CALL(sleqp_sparse_vector_create_empty(&stationarity_residuals,
                                               constrained_num_variables));

  ASSERT_CALL(sleqp_solver_solve(solver, 100, -1));

  ASSERT_CALL(
    sleqp_solver_get_vec_state(solver,
                               SLEQP_SOLVER_STATE_VEC_SCALED_STAT_RESIDUALS,
                               stationarity_residuals));

  const double stat_eps
    = sleqp_params_get(params, SLEQP_PARAM_STATIONARITY_TOL);

  const double eps = sleqp_params_get(params, SLEQP_PARAM_EPS);

  ck_assert(sleqp_is_leq(sleqp_sparse_vector_inf_norm(stationarity_residuals),
                         stat_eps,
                         eps));

  ASSERT_CALL(sleqp_sparse_vector_free(&stationarity_residuals));
}
END_TEST

START_TEST(test_feasibility_residuals)
{
  SleqpSparseVec* feasibility_residuals;

  ASSERT_CALL(sleqp_sparse_vector_create_empty(&feasibility_residuals,
                                               constrained_num_constraints));

  ASSERT_CALL(sleqp_solver_solve(solver, 100, -1));

  ASSERT_CALL(
    sleqp_solver_get_vec_state(solver,
                               SLEQP_SOLVER_STATE_VEC_SCALED_FEAS_RESIDUALS,
                               feasibility_residuals));

  const double feas_eps = sleqp_params_get(params, SLEQP_PARAM_FEASIBILITY_TOL);

  const double eps = sleqp_params_get(params, SLEQP_PARAM_EPS);

  ck_assert(sleqp_is_leq(sleqp_sparse_vector_inf_norm(feasibility_residuals),
                         feas_eps,
                         eps));

  ASSERT_CALL(sleqp_sparse_vector_free(&feasibility_residuals));
}
END_TEST

START_TEST(test_slackness_residuals)
{
  SleqpSparseVec* slackness_residuals;

  ASSERT_CALL(sleqp_sparse_vector_create_empty(&slackness_residuals,
                                               constrained_num_constraints));

  ASSERT_CALL(sleqp_solver_solve(solver, 100, -1));

  ASSERT_CALL(sleqp_solver_get_vec_state(
    solver,
    SLEQP_SOLVER_STATE_VEC_SCALED_CONS_SLACK_RESIDUALS,
    slackness_residuals));

  const double slack_eps = sleqp_params_get(params, SLEQP_PARAM_SLACKNESS_TOL);

  const double eps = sleqp_params_get(params, SLEQP_PARAM_EPS);

  ck_assert(sleqp_is_leq(sleqp_sparse_vector_inf_norm(slackness_residuals),
                         slack_eps,
                         eps));

  ASSERT_CALL(sleqp_sparse_vector_free(&slackness_residuals));
}
END_TEST

START_TEST(test_trust_radii)
{
  ASSERT_CALL(sleqp_solver_solve(solver, 100, -1));

  double trust_radius;
  double lp_trust_radius;

  ASSERT_CALL(sleqp_solver_get_real_state(solver,
                                          SLEQP_SOLVER_STATE_REAL_TRUST_RADIUS,
                                          &trust_radius));

  ck_assert(trust_radius > 0.);

  ASSERT_CALL(
    sleqp_solver_get_real_state(solver,
                                SLEQP_SOLVER_STATE_REAL_LP_TRUST_RADIUS,
                                &lp_trust_radius));

  ck_assert(lp_trust_radius > 0.);
}
END_TEST

START_TEST(test_residuals)
{
  ASSERT_CALL(sleqp_solver_solve(solver, 100, -1));

  const double eps = sleqp_params_get(params, SLEQP_PARAM_EPS);

  const double stat_eps
    = sleqp_params_get(params, SLEQP_PARAM_STATIONARITY_TOL);
  const double slack_eps = sleqp_params_get(params, SLEQP_PARAM_SLACKNESS_TOL);
  const double feas_eps = sleqp_params_get(params, SLEQP_PARAM_FEASIBILITY_TOL);

  double stat_res, slack_res, feas_res;

  ASSERT_CALL(
    sleqp_solver_get_real_state(solver,
                                SLEQP_SOLVER_STATE_REAL_SCALED_STAT_RES,
                                &stat_res));

  ck_assert(sleqp_is_leq(stat_res, stat_eps, eps));

  ASSERT_CALL(
    sleqp_solver_get_real_state(solver,
                                SLEQP_SOLVER_STATE_REAL_SCALED_FEAS_RES,
                                &feas_res));

  ck_assert(sleqp_is_leq(feas_res, feas_eps, eps));

  ASSERT_CALL(
    sleqp_solver_get_real_state(solver,
                                SLEQP_SOLVER_STATE_REAL_SCALED_SLACK_RES,
                                &slack_res));

  ck_assert(sleqp_is_leq(slack_res, slack_eps, eps));
}
END_TEST

START_TEST(test_iteration)
{
  ASSERT_CALL(sleqp_solver_solve(solver, 100, -1));

  int iteration;

  ASSERT_CALL(sleqp_solver_get_int_state(solver,
                                         SLEQP_SOLVER_STATE_INT_ITERATION,
                                         &iteration));

  ck_assert_int_gt(iteration, 0);
}
END_TEST

START_TEST(test_step_type)
{
  int last_step;

  ASSERT_CALL(sleqp_solver_get_int_state(solver,
                                         SLEQP_SOLVER_STATE_INT_LAST_STEP_TYPE,
                                         &last_step));

  ck_assert_int_ge(last_step, SLEQP_STEPTYPE_NONE);

  ck_assert_int_le(last_step, SLEQP_STEPTYPE_REJECTED);
}
END_TEST

Suite*
solver_state_test_suite()
{
  Suite* suite;
  TCase* tc_state;

  suite = suite_create("Solver state tests");

  tc_state = tcase_create("Solver state test");

  tcase_add_checked_fixture(tc_state,
                            solver_state_setup,
                            solver_state_teardown);

  tcase_add_test(tc_state, test_stationarity_residuals);

  tcase_add_test(tc_state, test_feasibility_residuals);

  tcase_add_test(tc_state, test_slackness_residuals);

  tcase_add_test(tc_state, test_trust_radii);

  tcase_add_test(tc_state, test_residuals);

  tcase_add_test(tc_state, test_iteration);

  tcase_add_test(tc_state, test_step_type);

  suite_add_tcase(suite, tc_state);

  return suite;
}

TEST_MAIN(solver_state_test_suite)
