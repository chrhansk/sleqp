#include <check.h>
#include <stdlib.h>

#include "cmp.h"
#include "mem.h"
#include "solver.h"
#include "util.h"

#include "constrained_fixture.h"
#include "test_common.h"

SleqpSettings* settings;
SleqpProblem* problem;

void
constrained_test_setup()
{
  constrained_setup();

  ASSERT_CALL(sleqp_settings_create(&settings));

  ASSERT_CALL(sleqp_problem_create_simple(&problem,
                                          constrained_func,
                                          constrained_var_lb,
                                          constrained_var_ub,
                                          constrained_cons_lb,
                                          constrained_cons_ub,
                                          settings));
}

void
constrained_test_teardown()
{
  ASSERT_CALL(sleqp_problem_release(&problem));

  ASSERT_CALL(sleqp_settings_release(&settings));

  constrained_teardown();
}

double
stat_res(SleqpProblem* problem, SleqpIterate* iterate)
{
  const int num_vars = sleqp_problem_num_vars(problem);

  double* residuals;
  SleqpVec* product;

  ASSERT_CALL(sleqp_alloc_array(&residuals, num_vars));
  ASSERT_CALL(sleqp_vec_create_full(&product, num_vars));

  double stat_res = 0.;

  ASSERT_CALL(sleqp_vec_to_raw(sleqp_iterate_obj_grad(iterate), residuals));

  ASSERT_CALL(sleqp_mat_mult_vec_trans(sleqp_iterate_cons_jac(iterate),
                                       sleqp_iterate_cons_dual(iterate),
                                       0.,
                                       product));

  for (int k = 0; k < product->nnz; ++k)
  {
    residuals[product->indices[k]] += product->data[k];
  }

  SleqpVec* vars_dual = sleqp_iterate_vars_dual(iterate);

  for (int k = 0; k < vars_dual->nnz; ++k)
  {
    residuals[vars_dual->indices[k]] += vars_dual->data[k];
  }

  for (int i = 0; i < num_vars; ++i)
  {
    stat_res = SLEQP_MAX(stat_res, residuals[i]);
  }

  ASSERT_CALL(sleqp_vec_free(&product));
  sleqp_free(&residuals);

  return stat_res;
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

  ck_assert(sleqp_is_zero(stat_res(problem, iterate), 1e-6));

  ASSERT_CALL(sleqp_solver_release(&solver));
}

START_TEST(test_solve)
{
  SleqpSolver* solver;

  ASSERT_CALL(sleqp_settings_set_enum_value(
    settings,
    SLEQP_SETTINGS_ENUM_DERIV_CHECK,
    SLEQP_DERIV_CHECK_FIRST | SLEQP_DERIV_CHECK_SECOND_EXHAUSTIVE));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  constrained_initial,
                                  NULL));

  solve_and_release_solver(solver);
}
END_TEST

START_TEST(test_solve_reduced)
{
  SleqpSolver* solver;

  ASSERT_CALL(sleqp_settings_set_enum_value(settings,
                                           SLEQP_SETTINGS_ENUM_AUG_JAC_METHOD,
                                           SLEQP_AUG_JAC_REDUCED));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  constrained_initial,
                                  NULL));

  solve_and_release_solver(solver);
}
END_TEST

#ifdef SLEQP_HAVE_QR_FACT

START_TEST(test_solve_direct)
{
  SleqpSolver* solver;

  ASSERT_CALL(sleqp_settings_set_enum_value(settings,
                                           SLEQP_SETTINGS_ENUM_AUG_JAC_METHOD,
                                           SLEQP_AUG_JAC_DIRECT));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  constrained_initial,
                                  NULL));

  solve_and_release_solver(solver);
}
END_TEST

#endif

START_TEST(test_exact_linesearch)
{
  SleqpSolver* solver;

  ASSERT_CALL(sleqp_settings_set_enum_value(settings,
                                           SLEQP_SETTINGS_ENUM_LINESEARCH,
                                           SLEQP_LINESEARCH_EXACT));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  constrained_initial,
                                  NULL));

  solve_and_release_solver(solver);
}
END_TEST

START_TEST(test_initial_tr_wide)
{
  SleqpSolver* solver;

  ASSERT_CALL(sleqp_settings_set_enum_value(settings,
                                           SLEQP_SETTINGS_ENUM_INITIAL_TR_CHOICE,
                                           SLEQP_INITIAL_TR_CHOICE_WIDE));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  constrained_initial,
                                  NULL));

  solve_and_release_solver(solver);
}
END_TEST

START_TEST(test_parametric_solve)
{
  SleqpSolver* solver;

  ASSERT_CALL(sleqp_settings_set_enum_value(settings,
                                           SLEQP_SETTINGS_ENUM_PARAMETRIC_CAUCHY,
                                           SLEQP_PARAMETRIC_CAUCHY_COARSE));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  constrained_initial,
                                  NULL));

  solve_and_release_solver(solver);
}
END_TEST

START_TEST(test_sr1_solve)
{
  ASSERT_CALL(sleqp_settings_set_enum_value(settings,
                                           SLEQP_SETTINGS_ENUM_DERIV_CHECK,
                                           SLEQP_DERIV_CHECK_FIRST));

  ASSERT_CALL(sleqp_settings_set_enum_value(settings,
                                           SLEQP_SETTINGS_ENUM_HESS_EVAL,
                                           SLEQP_HESS_EVAL_SR1));

  SleqpSolver* solver;

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  constrained_initial,
                                  NULL));

  solve_and_release_solver(solver);
}
END_TEST

START_TEST(test_bfgs_solve_no_sizing)
{
  ASSERT_CALL(sleqp_settings_set_enum_value(settings,
                                           SLEQP_SETTINGS_ENUM_DERIV_CHECK,
                                           SLEQP_DERIV_CHECK_FIRST));

  ASSERT_CALL(sleqp_settings_set_enum_value(settings,
                                           SLEQP_SETTINGS_ENUM_BFGS_SIZING,
                                           SLEQP_BFGS_SIZING_NONE));

  ASSERT_CALL(sleqp_settings_set_enum_value(settings,
                                           SLEQP_SETTINGS_ENUM_HESS_EVAL,
                                           SLEQP_HESS_EVAL_DAMPED_BFGS));

  SleqpSolver* solver;

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  constrained_initial,
                                  NULL));

  solve_and_release_solver(solver);
}
END_TEST

START_TEST(test_bfgs_solve_centered_ol_sizing)
{
  ASSERT_CALL(sleqp_settings_set_enum_value(settings,
                                           SLEQP_SETTINGS_ENUM_BFGS_SIZING,
                                           SLEQP_BFGS_SIZING_CENTERED_OL));

  ASSERT_CALL(sleqp_settings_set_enum_value(settings,
                                           SLEQP_SETTINGS_ENUM_HESS_EVAL,
                                           SLEQP_HESS_EVAL_DAMPED_BFGS));

  SleqpSolver* solver;

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  constrained_initial,
                                  NULL));

  solve_and_release_solver(solver);
}
END_TEST

START_TEST(test_unscaled_solve)
{
  ASSERT_CALL(sleqp_settings_set_enum_value(settings,
                                           SLEQP_SETTINGS_ENUM_DERIV_CHECK,
                                           SLEQP_DERIV_CHECK_FIRST));

  SleqpSolver* solver;

  SleqpScaling* scaling;

  ASSERT_CALL(sleqp_scaling_create(&scaling,
                                   constrained_num_variables,
                                   constrained_num_constraints));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  constrained_initial,
                                  scaling));

  solve_and_release_solver(solver);

  ASSERT_CALL(sleqp_scaling_release(&scaling));
}
END_TEST

START_TEST(test_scaled_solve)
{
  ASSERT_CALL(sleqp_settings_set_enum_value(
    settings,
    SLEQP_SETTINGS_ENUM_DERIV_CHECK,
    SLEQP_DERIV_CHECK_FIRST | SLEQP_DERIV_CHECK_SECOND_EXHAUSTIVE));

  ASSERT_CALL(sleqp_settings_set_real_value(settings, SLEQP_SETTINGS_REAL_STAT_TOL, 1e-7));

  SleqpSolver* solver;

  SleqpScaling* scaling;

  ASSERT_CALL(sleqp_scaling_create(&scaling,
                                   constrained_num_variables,
                                   constrained_num_constraints));

  ASSERT_CALL(sleqp_scaling_set_obj_weight(scaling, 2));

  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling, 0, -5));
  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling, 1, 5));

  ASSERT_CALL(sleqp_scaling_set_cons_weight(scaling, 0, -10));
  ASSERT_CALL(sleqp_scaling_set_cons_weight(scaling, 1, 10));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  constrained_initial,
                                  scaling));

  solve_and_release_solver(solver);

  ASSERT_CALL(sleqp_scaling_release(&scaling));
}
END_TEST

START_TEST(test_scaled_sr1_solve)
{
  ASSERT_CALL(sleqp_settings_set_enum_value(settings,
                                           SLEQP_SETTINGS_ENUM_DERIV_CHECK,
                                           SLEQP_DERIV_CHECK_FIRST));

  ASSERT_CALL(sleqp_settings_set_real_value(settings, SLEQP_SETTINGS_REAL_STAT_TOL, 1e-7));

  SleqpSolver* solver;

  SleqpScaling* scaling;

  ASSERT_CALL(sleqp_settings_set_enum_value(settings,
                                           SLEQP_SETTINGS_ENUM_HESS_EVAL,
                                           SLEQP_HESS_EVAL_SR1));

  ASSERT_CALL(sleqp_scaling_create(&scaling,
                                   constrained_num_variables,
                                   constrained_num_constraints));

  ASSERT_CALL(sleqp_scaling_set_obj_weight(scaling, 2));

  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling, 0, -5));
  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling, 1, 5));

  ASSERT_CALL(sleqp_scaling_set_cons_weight(scaling, 0, -1));
  ASSERT_CALL(sleqp_scaling_set_cons_weight(scaling, 1, 2));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  constrained_initial,
                                  scaling));

  solve_and_release_solver(solver);

  ASSERT_CALL(sleqp_scaling_release(&scaling));
}
END_TEST

START_TEST(test_scaled_bfgs_solve)
{
  ASSERT_CALL(sleqp_settings_set_enum_value(settings,
                                           SLEQP_SETTINGS_ENUM_DERIV_CHECK,
                                           SLEQP_DERIV_CHECK_FIRST));

  ASSERT_CALL(sleqp_settings_set_real_value(settings, SLEQP_SETTINGS_REAL_STAT_TOL, 1e-7));

  SleqpSolver* solver;

  SleqpScaling* scaling;

  ASSERT_CALL(sleqp_settings_set_enum_value(settings,
                                           SLEQP_SETTINGS_ENUM_HESS_EVAL,
                                           SLEQP_HESS_EVAL_DAMPED_BFGS));

  ASSERT_CALL(sleqp_scaling_create(&scaling,
                                   constrained_num_variables,
                                   constrained_num_constraints));

  ASSERT_CALL(sleqp_scaling_set_obj_weight(scaling, 2));

  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling, 0, -5));
  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling, 1, 5));

  ASSERT_CALL(sleqp_scaling_set_cons_weight(scaling, 0, -1));
  ASSERT_CALL(sleqp_scaling_set_cons_weight(scaling, 1, 2));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  constrained_initial,
                                  scaling));

  solve_and_release_solver(solver);

  ASSERT_CALL(sleqp_scaling_release(&scaling));
}
END_TEST

START_TEST(test_auto_scaled_solve)
{
  ASSERT_CALL(sleqp_settings_set_enum_value(settings,
                                           SLEQP_SETTINGS_ENUM_DERIV_CHECK,
                                           SLEQP_DERIV_CHECK_FIRST));

  SleqpSolver* solver;

  SleqpScaling* scaling;

  SleqpIterate* iterate;

  const double eps = sleqp_settings_real_value(settings, SLEQP_SETTINGS_REAL_EPS);

  ASSERT_CALL(sleqp_iterate_create(&iterate, problem, constrained_initial));

  ASSERT_CALL(
    sleqp_set_and_evaluate(problem, iterate, SLEQP_VALUE_REASON_NONE, NULL));

  ASSERT_CALL(sleqp_scaling_create(&scaling,
                                   constrained_num_variables,
                                   constrained_num_constraints));

  ASSERT_CALL(
    sleqp_obj_scaling_from_grad(scaling, sleqp_iterate_obj_grad(iterate), eps));

  ASSERT_CALL(
    sleqp_scaling_from_cons_jac(scaling, sleqp_iterate_cons_jac(iterate), eps));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  constrained_initial,
                                  scaling));

  solve_and_release_solver(solver);

  ASSERT_CALL(sleqp_scaling_release(&scaling));

  ASSERT_CALL(sleqp_iterate_release(&iterate));
}
END_TEST

START_TEST(test_lp_dual_estimation)
{
  SleqpSolver* solver;

  ASSERT_CALL(
    sleqp_settings_set_enum_value(settings,
                                 SLEQP_SETTINGS_ENUM_DUAL_ESTIMATION_TYPE,
                                 SLEQP_DUAL_ESTIMATION_TYPE_LP));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  constrained_initial,
                                  NULL));

  solve_and_release_solver(solver);
}
END_TEST

START_TEST(test_mixed_dual_estimation)
{
  SleqpSolver* solver;

  ASSERT_CALL(
    sleqp_settings_set_enum_value(settings,
                                 SLEQP_SETTINGS_ENUM_DUAL_ESTIMATION_TYPE,
                                 SLEQP_DUAL_ESTIMATION_TYPE_MIXED));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  constrained_initial,
                                  NULL));

  solve_and_release_solver(solver);
}
END_TEST

Suite*
constrained_test_suite()
{
  Suite* suite;
  TCase* tc_cons;

  suite = suite_create("Constrained tests");

  tc_cons = tcase_create("Constrained solution test");

  tcase_add_checked_fixture(tc_cons,
                            constrained_test_setup,
                            constrained_test_teardown);

  tcase_add_test(tc_cons, test_solve);

  tcase_add_test(tc_cons, test_solve_reduced);

#ifdef SLEQP_HAVE_QR_FACT

  tcase_add_test(tc_cons, test_solve_direct);

#endif

  tcase_add_test(tc_cons, test_exact_linesearch);

  tcase_add_test(tc_cons, test_initial_tr_wide);

  tcase_add_test(tc_cons, test_parametric_solve);

  tcase_add_test(tc_cons, test_sr1_solve);

  tcase_add_test(tc_cons, test_bfgs_solve_no_sizing);

  tcase_add_test(tc_cons, test_bfgs_solve_centered_ol_sizing);

  tcase_add_test(tc_cons, test_unscaled_solve);

  tcase_add_test(tc_cons, test_scaled_solve);

  tcase_add_test(tc_cons, test_scaled_sr1_solve);

  tcase_add_test(tc_cons, test_scaled_bfgs_solve);

  tcase_add_test(tc_cons, test_auto_scaled_solve);

  tcase_add_test(tc_cons, test_lp_dual_estimation);

  tcase_add_test(tc_cons, test_mixed_dual_estimation);

  suite_add_tcase(suite, tc_cons);

  return suite;
}

TEST_MAIN(constrained_test_suite)
