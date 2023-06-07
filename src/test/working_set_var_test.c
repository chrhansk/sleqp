#include <check.h>
#include <stdlib.h>

#include "cauchy/standard_cauchy.h"
#include "cmp.h"
#include "lp/lpi.h"
#include "mem.h"
#include "util.h"

#include "quadfunc_fixture.h"
#include "test_common.h"

SleqpSettings* settings;
SleqpProblem* problem;
SleqpIterate* iterate;
SleqpCauchy* cauchy_data;

void
working_set_var_setup()
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

  ASSERT_CALL(
    sleqp_set_and_evaluate(problem, iterate, SLEQP_VALUE_REASON_NONE, NULL));

  ASSERT_CALL(
    sleqp_standard_cauchy_create(&cauchy_data, problem, settings));
}

START_TEST(test_inactive)
{
  SleqpVec* primal = sleqp_iterate_primal(iterate);

  primal->data[0] = 1.5;
  primal->data[1] = 2.5;

  ASSERT_CALL(
    sleqp_set_and_evaluate(problem, iterate, SLEQP_VALUE_REASON_NONE, NULL));

  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);

  double trust_radius = 0.25, penalty_parameter = 1.;

  ASSERT_CALL(sleqp_cauchy_set_iterate(cauchy_data, iterate, trust_radius));

  ASSERT_CALL(sleqp_cauchy_solve(cauchy_data,
                                 sleqp_iterate_obj_grad(iterate),
                                 penalty_parameter,
                                 SLEQP_CAUCHY_OBJTYPE_DEFAULT));

  ASSERT_CALL(sleqp_cauchy_working_set(cauchy_data, iterate));

  ck_assert_int_eq(sleqp_working_set_var_state(working_set, 0), SLEQP_INACTIVE);
  ck_assert_int_eq(sleqp_working_set_var_state(working_set, 1), SLEQP_INACTIVE);
}
END_TEST

START_TEST(test_active)
{
  SleqpVec* primal = sleqp_iterate_primal(iterate);

  primal->data[0] = 1.5;
  primal->data[1] = 2.5;

  ASSERT_CALL(
    sleqp_set_and_evaluate(problem, iterate, SLEQP_VALUE_REASON_NONE, NULL));

  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);

  double trust_radius = 1., penalty_parameter = 1.;

  ASSERT_CALL(sleqp_cauchy_set_iterate(cauchy_data, iterate, trust_radius));

  ASSERT_CALL(sleqp_cauchy_solve(cauchy_data,
                                 sleqp_iterate_obj_grad(iterate),
                                 penalty_parameter,
                                 SLEQP_CAUCHY_OBJTYPE_DEFAULT));

  ASSERT_CALL(sleqp_cauchy_working_set(cauchy_data, iterate));

  ck_assert_int_eq(sleqp_working_set_var_state(working_set, 0),
                   SLEQP_ACTIVE_LOWER);
  ck_assert_int_eq(sleqp_working_set_var_state(working_set, 1),
                   SLEQP_ACTIVE_LOWER);
}
END_TEST

START_TEST(test_first_active)
{
  SleqpVec* primal = sleqp_iterate_primal(iterate);

  primal->data[0] = 1.;
  primal->data[1] = 2.5;

  ASSERT_CALL(
    sleqp_set_and_evaluate(problem, iterate, SLEQP_VALUE_REASON_NONE, NULL));

  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);

  double trust_radius = 0.25, penalty_parameter = 1.;

  ASSERT_CALL(sleqp_cauchy_set_iterate(cauchy_data, iterate, trust_radius));

  ASSERT_CALL(sleqp_cauchy_solve(cauchy_data,
                                 sleqp_iterate_obj_grad(iterate),
                                 penalty_parameter,
                                 SLEQP_CAUCHY_OBJTYPE_DEFAULT));

  ASSERT_CALL(sleqp_cauchy_working_set(cauchy_data, iterate));

  ck_assert_int_eq(sleqp_working_set_var_state(working_set, 0),
                   SLEQP_ACTIVE_LOWER);
  ck_assert_int_eq(sleqp_working_set_var_state(working_set, 1), SLEQP_INACTIVE);
}
END_TEST

void
working_set_var_teardown()
{
  ASSERT_CALL(sleqp_cauchy_release(&cauchy_data));

  ASSERT_CALL(sleqp_iterate_release(&iterate));

  ASSERT_CALL(sleqp_problem_release(&problem));

  ASSERT_CALL(sleqp_settings_release(&settings));

  quadfunc_teardown();
}

Suite*
working_set_var_test_suite()
{
  Suite* suite;
  TCase* tc_working_set_var;

  suite = suite_create("Dual estimation tests");

  tc_working_set_var = tcase_create("Simply constrained");

  tcase_add_checked_fixture(tc_working_set_var,
                            working_set_var_setup,
                            working_set_var_teardown);

  tcase_add_test(tc_working_set_var, test_inactive);
  tcase_add_test(tc_working_set_var, test_active);
  tcase_add_test(tc_working_set_var, test_first_active);

  suite_add_tcase(suite, tc_working_set_var);

  return suite;
}

TEST_MAIN(working_set_var_test_suite)
