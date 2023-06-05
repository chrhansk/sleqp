#include <check.h>
#include <stdlib.h>

#include "aug_jac/standard_aug_jac.h"
#include "cauchy/standard_cauchy.h"
#include "cmp.h"
#include "direction.h"
#include "dual_estimation/dual_estimation.h"
#include "fact/fact.h"
#include "mem.h"
#include "newton.h"
#include "pub_types.h"
#include "util.h"
#include "working_step.h"

#include "quadfunc_fixture.h"
#include "test_common.h"

SleqpSettings* settings;
SleqpProblem* problem;
SleqpIterate* iterate;
SleqpCauchy* cauchy_data;

const double tolerance = 1e-8;

void
newton_setup()
{
  quadfunc_setup();

  ASSERT_CALL(sleqp_settings_create(&settings));

  ASSERT_CALL(sleqp_problem_create_simple(&problem,
                                          quadfunc,
                                          quadfunc_var_lb,
                                          quadfunc_var_ub,
                                          quadfunc_cons_lb,
                                          quadfunc_cons_ub,
                                          settings));

  ASSERT_CALL(sleqp_iterate_create(&iterate, problem, quadfunc_x));

  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);
  ASSERT_CALL(sleqp_working_set_reset(working_set));

  ASSERT_CALL(
    sleqp_set_and_evaluate(problem, iterate, SLEQP_VALUE_REASON_INIT, NULL));

  ASSERT_CALL(
    sleqp_standard_cauchy_create(&cauchy_data, problem, settings));
}

void
newton_teardown()
{
  ASSERT_CALL(sleqp_cauchy_release(&cauchy_data));

  ASSERT_CALL(sleqp_iterate_release(&iterate));

  ASSERT_CALL(sleqp_problem_release(&problem));

  ASSERT_CALL(sleqp_settings_release(&settings));

  quadfunc_teardown();
}

START_TEST(newton_wide_step)
{
  SleqpWorkingStep* working_step;
  SleqpEQPSolver* newton_solver;

  SleqpVec* expected_step;
  SleqpDirection* actual_direction;

  SleqpFact* fact;
  SleqpAugJac* jacobian;

  const int num_variables = sleqp_problem_num_vars(problem);

  ASSERT_CALL(sleqp_vec_create(&expected_step, num_variables, 2));

  ASSERT_CALL(sleqp_vec_push(expected_step, 0, -1.));
  ASSERT_CALL(sleqp_vec_push(expected_step, 1, -2.));

  // ASSERT_CALL(sleqp_vec_create(&actual_step, num_variables, 0));
  ASSERT_CALL(sleqp_direction_create(&actual_direction, problem, settings));

  double penalty_parameter = 1.;
  double trust_radius      = 10.;

  ASSERT_CALL(sleqp_fact_create_default(&fact, settings));

  // create with empty active set
  ASSERT_CALL(sleqp_standard_aug_jac_create(&jacobian, problem, settings, fact));

  ASSERT_CALL(sleqp_aug_jac_set_iterate(jacobian, iterate));

  ASSERT_CALL(sleqp_working_step_create(&working_step, problem, settings));

  ASSERT_CALL(sleqp_newton_solver_create(&newton_solver,
                                         problem,
                                         settings,
                                         working_step));

  ASSERT_CALL(sleqp_eqp_solver_set_iterate(newton_solver,
                                           iterate,
                                           jacobian,
                                           trust_radius,
                                           penalty_parameter));

  // we use the default (empty) active set for the Newton step,
  // trust region size should be large to ensure that
  // the solution is that of the unrestricted step
  ASSERT_CALL(
    sleqp_eqp_solver_compute_direction(newton_solver,
                                       sleqp_iterate_cons_dual(iterate),
                                       actual_direction));

  SleqpVec* actual_step = sleqp_direction_primal(actual_direction);

  ck_assert(sleqp_vec_eq(expected_step, actual_step, tolerance));

  ASSERT_CALL(sleqp_eqp_solver_release(&newton_solver));

  ASSERT_CALL(sleqp_working_step_release(&working_step));

  ASSERT_CALL(sleqp_aug_jac_release(&jacobian));

  ASSERT_CALL(sleqp_fact_release(&fact));

  ASSERT_CALL(sleqp_direction_release(&actual_direction));

  ASSERT_CALL(sleqp_vec_free(&expected_step));
}
END_TEST

START_TEST(newton_small_step)
{
  SleqpWorkingStep* working_step;
  SleqpEQPSolver* newton_solver;

  SleqpVec* expected_step;
  SleqpDirection* actual_direction;

  ASSERT_CALL(sleqp_direction_create(&actual_direction, problem, settings));

  SleqpFact* factorization;
  SleqpAugJac* jacobian;

  const int num_variables = sleqp_problem_num_vars(problem);

  ASSERT_CALL(sleqp_vec_create(&expected_step, num_variables, 2));

  ASSERT_CALL(sleqp_vec_push(expected_step, 0, -0.44721359549995793));
  ASSERT_CALL(sleqp_vec_push(expected_step, 1, -0.89442719099991586));

  double penalty_parameter = 1.;
  double trust_radius      = 1.;

  ASSERT_CALL(sleqp_fact_create_default(&factorization, settings));

  // create with empty active set
  ASSERT_CALL(
    sleqp_standard_aug_jac_create(&jacobian, problem, settings, factorization));

  ASSERT_CALL(sleqp_aug_jac_set_iterate(jacobian, iterate));

  ASSERT_CALL(sleqp_working_step_create(&working_step, problem, settings));

  ASSERT_CALL(sleqp_newton_solver_create(&newton_solver,
                                         problem,
                                         settings,
                                         working_step));

  ASSERT_CALL(sleqp_eqp_solver_set_iterate(newton_solver,
                                           iterate,
                                           jacobian,
                                           trust_radius,
                                           penalty_parameter));

  // we use the default (empty) active set for the Newton step,
  // trust region size should be so small that
  // the solution is on the boundary of the feasible set
  ASSERT_CALL(
    sleqp_eqp_solver_compute_direction(newton_solver,
                                       sleqp_iterate_cons_dual(iterate),
                                       actual_direction));

  SleqpVec* actual_step = sleqp_direction_primal(actual_direction);

  ck_assert(sleqp_vec_eq(expected_step, actual_step, tolerance));

  ASSERT_CALL(sleqp_eqp_solver_release(&newton_solver));

  ASSERT_CALL(sleqp_working_step_release(&working_step));

  ASSERT_CALL(sleqp_aug_jac_release(&jacobian));

  ASSERT_CALL(sleqp_fact_release(&factorization));

  ASSERT_CALL(sleqp_direction_release(&actual_direction));

  ASSERT_CALL(sleqp_vec_free(&expected_step));
}
END_TEST

Suite*
newton_test_suite()
{
  Suite* suite;
  TCase* tc_cons;

  suite = suite_create("Unconstrained newton step tests");

  tc_cons = tcase_create("Newton step");

  tcase_add_checked_fixture(tc_cons, newton_setup, newton_teardown);

  tcase_add_test(tc_cons, newton_wide_step);

  tcase_add_test(tc_cons, newton_small_step);

  suite_add_tcase(suite, tc_cons);

  return suite;
}

TEST_MAIN(newton_test_suite)
