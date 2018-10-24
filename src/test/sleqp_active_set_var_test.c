#include <stdlib.h>
#include <check.h>

#include "sleqp.h"
#include "sleqp_aug_jacobian.h"
#include "sleqp_cauchy.h"
#include "sleqp_cmp.h"
#include "sleqp_dual_estimation.h"
#include "sleqp_mem.h"

#include "lp/sleqp_lpi_soplex.h"

#include "test_common.h"

#include "quadfunc_fixture.h"

SleqpProblem* problem;
SleqpIterate* iterate;
SleqpLPi* lp_interface;
SleqpCauchyData* cauchy_data;

void active_set_var_setup()
{
  quadfunc_setup();

  ASSERT_CALL(sleqp_problem_create(&problem,
                                   quadfunc,
                                   quadfunc_var_lb,
                                   quadfunc_var_ub,
                                   quadfunc_cons_lb,
                                   quadfunc_cons_ub));

  ASSERT_CALL(sleqp_iterate_create(&iterate,
                                   problem,
                                   quadfunc_x));

  int num_lp_variables = problem->num_variables + 2*problem->num_constraints;
  int num_lp_constraints = problem->num_constraints;

  ASSERT_CALL(sleqp_lpi_soplex_create_interface(&lp_interface,
                                                num_lp_variables,
                                                num_lp_constraints));

  ASSERT_CALL(sleqp_set_and_evaluate(problem, iterate));

  ASSERT_CALL(sleqp_cauchy_data_create(&cauchy_data,
                                       problem,
                                       lp_interface));
}

START_TEST(test_inactive)
{
  iterate->x->data[0] = 1.5;
  iterate->x->data[1] = 2.5;

  ASSERT_CALL(sleqp_set_and_evaluate(problem, iterate));

  SleqpActiveSet* active_set = iterate->active_set;

  double trust_radius = 0.25, penalty = 1.;

  ASSERT_CALL(sleqp_cauchy_compute_direction(cauchy_data,
                                             iterate,
                                             penalty,
                                             trust_radius));

  ASSERT_CALL(sleqp_cauchy_get_active_set(cauchy_data,
                                          iterate,
                                          trust_radius));

  SLEQP_ACTIVE_STATE* var_states = sleqp_active_set_var_states(active_set);

  ck_assert_int_eq(var_states[0], SLEQP_INACTIVE);
  ck_assert_int_eq(var_states[1], SLEQP_INACTIVE);
}
END_TEST

START_TEST(test_active)
{
  iterate->x->data[0] = 1.5;
  iterate->x->data[1] = 2.5;

  ASSERT_CALL(sleqp_set_and_evaluate(problem, iterate));

  SleqpActiveSet* active_set = iterate->active_set;

  double trust_radius = 1., penalty = 1.;

  ASSERT_CALL(sleqp_cauchy_compute_direction(cauchy_data,
                                             iterate,
                                             penalty,
                                             trust_radius));

  ASSERT_CALL(sleqp_cauchy_get_active_set(cauchy_data,
                                          iterate,
                                          trust_radius));

  SLEQP_ACTIVE_STATE* var_states = sleqp_active_set_var_states(active_set);

  ck_assert_int_eq(var_states[0], SLEQP_ACTIVE_LOWER);
  ck_assert_int_eq(var_states[1], SLEQP_ACTIVE_LOWER);
}
END_TEST

START_TEST(test_first_active)
{
  iterate->x->data[0] = 1.;
  iterate->x->data[1] = 2.5;

  ASSERT_CALL(sleqp_set_and_evaluate(problem, iterate));

  SleqpActiveSet* active_set = iterate->active_set;

  double trust_radius = 0.25, penalty = 1.;

  ASSERT_CALL(sleqp_cauchy_compute_direction(cauchy_data,
                                             iterate,
                                             penalty,
                                             trust_radius));

  ASSERT_CALL(sleqp_cauchy_get_active_set(cauchy_data,
                                          iterate,
                                          trust_radius));

  SLEQP_ACTIVE_STATE* var_states = sleqp_active_set_var_states(active_set);

  ck_assert_int_eq(var_states[0], SLEQP_ACTIVE_LOWER);
  ck_assert_int_eq(var_states[1], SLEQP_INACTIVE);
}
END_TEST

void active_set_var_teardown()
{
  ASSERT_CALL(sleqp_cauchy_data_free(&cauchy_data));

  ASSERT_CALL(sleqp_lpi_free(&lp_interface));

  ASSERT_CALL(sleqp_iterate_free(&iterate));

  ASSERT_CALL(sleqp_problem_free(&problem));

  quadfunc_teardown();
}

Suite* active_set_var_test_suite()
{
  Suite *suite;
  TCase *tc_active_set_var;

  suite = suite_create("Dual estimation tests");

  tc_active_set_var = tcase_create("Simply constrained");

  tcase_add_checked_fixture(tc_active_set_var,
                            active_set_var_setup,
                            active_set_var_teardown);

  tcase_add_test(tc_active_set_var, test_inactive);
  tcase_add_test(tc_active_set_var, test_active);
  tcase_add_test(tc_active_set_var, test_first_active);

  suite_add_tcase(suite, tc_active_set_var);

  return suite;
}


int main()
{
  int num_fails;
  Suite* suite;
  SRunner* srunner;

  suite = active_set_var_test_suite();
  srunner = srunner_create(suite);

  srunner_set_fork_status(srunner, CK_NOFORK);
  srunner_run_all(srunner, CK_NORMAL);

  num_fails = srunner_ntests_failed(srunner);

  srunner_free(srunner);

  return (num_fails > 0) ? EXIT_FAILURE : EXIT_SUCCESS;
}
