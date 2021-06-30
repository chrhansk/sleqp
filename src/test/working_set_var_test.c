#include <stdlib.h>
#include <check.h>

#include "sleqp.h"
#include "sleqp_aug_jacobian.h"
#include "sleqp_cauchy.h"
#include "sleqp_cmp.h"
#include "sleqp_dual_estimation.h"
#include "sleqp_mem.h"

#include "lp/sleqp_lpi.h"

#include "test_common.h"

#include "quadfunc_fixture.h"

SleqpParams* params;
SleqpOptions* options;
SleqpProblem* problem;
SleqpIterate* iterate;
SleqpLPi* lp_interface;
SleqpCauchy* cauchy_data;

void working_set_var_setup()
{
  quadfunc_setup();

  ASSERT_CALL(sleqp_params_create(&params));
  ASSERT_CALL(sleqp_options_create(&options));

  ASSERT_CALL(sleqp_problem_create_simple(&problem,
                                          quadfunc,
                                          params,
                                          quadfunc_var_lb,
                                          quadfunc_var_ub,
                                          quadfunc_cons_lb,
                                          quadfunc_cons_ub));

  const int num_variables = sleqp_problem_num_variables(problem);
  const int num_constraints = sleqp_problem_num_constraints(problem);

  ASSERT_CALL(sleqp_iterate_create(&iterate,
                                   problem,
                                   quadfunc_x));

  int num_lp_variables = num_variables + 2*num_constraints;
  int num_lp_constraints = num_constraints;

  ASSERT_CALL(sleqp_lpi_create_default_interface(&lp_interface,
                                                 num_lp_variables,
                                                 num_lp_constraints,
                                                 params));

  ASSERT_CALL(sleqp_set_and_evaluate(problem, iterate, SLEQP_VALUE_REASON_NONE));

  ASSERT_CALL(sleqp_cauchy_create(&cauchy_data,
                                  problem,
                                  params,
                                  options,
                                  lp_interface));
}

START_TEST(test_inactive)
{
  SleqpSparseVec* primal = sleqp_iterate_get_primal(iterate);

  primal->data[0] = 1.5;
  primal->data[1] = 2.5;

  ASSERT_CALL(sleqp_set_and_evaluate(problem, iterate, SLEQP_VALUE_REASON_NONE));

  SleqpWorkingSet* working_set = sleqp_iterate_get_working_set(iterate);

  double trust_radius = 0.25, penalty_parameter = 1.;

  ASSERT_CALL(sleqp_cauchy_set_iterate(cauchy_data,
                                       iterate,
                                       trust_radius));

  ASSERT_CALL(sleqp_cauchy_solve(cauchy_data,
                                 sleqp_iterate_get_func_grad(iterate),
                                 penalty_parameter,
                                 SLEQP_CAUCHY_OBJECTIVE_TYPE_DEFAULT));

  ASSERT_CALL(sleqp_cauchy_get_working_set(cauchy_data,
                                           iterate));

  ck_assert_int_eq(sleqp_working_set_get_variable_state(working_set, 0), SLEQP_INACTIVE);
  ck_assert_int_eq(sleqp_working_set_get_variable_state(working_set, 1), SLEQP_INACTIVE);
}
END_TEST

START_TEST(test_active)
{
  SleqpSparseVec* primal = sleqp_iterate_get_primal(iterate);

  primal->data[0] = 1.5;
  primal->data[1] = 2.5;

  ASSERT_CALL(sleqp_set_and_evaluate(problem, iterate, SLEQP_VALUE_REASON_NONE));

  SleqpWorkingSet* working_set = sleqp_iterate_get_working_set(iterate);

  double trust_radius = 1., penalty_parameter = 1.;

  ASSERT_CALL(sleqp_cauchy_set_iterate(cauchy_data,
                                       iterate,
                                       trust_radius));

  ASSERT_CALL(sleqp_cauchy_solve(cauchy_data,
                                 sleqp_iterate_get_func_grad(iterate),
                                 penalty_parameter,
                                 SLEQP_CAUCHY_OBJECTIVE_TYPE_DEFAULT));

  ASSERT_CALL(sleqp_cauchy_get_working_set(cauchy_data,
                                           iterate));

  ck_assert_int_eq(sleqp_working_set_get_variable_state(working_set, 0), SLEQP_ACTIVE_LOWER);
  ck_assert_int_eq(sleqp_working_set_get_variable_state(working_set, 1), SLEQP_ACTIVE_LOWER);
}
END_TEST

START_TEST(test_first_active)
{
  SleqpSparseVec* primal = sleqp_iterate_get_primal(iterate);

  primal->data[0] = 1.;
  primal->data[1] = 2.5;

  ASSERT_CALL(sleqp_set_and_evaluate(problem, iterate, SLEQP_VALUE_REASON_NONE));

  SleqpWorkingSet* working_set = sleqp_iterate_get_working_set(iterate);

  double trust_radius = 0.25, penalty_parameter = 1.;

  ASSERT_CALL(sleqp_cauchy_set_iterate(cauchy_data,
                                       iterate,
                                       trust_radius));

  ASSERT_CALL(sleqp_cauchy_solve(cauchy_data,
                                 sleqp_iterate_get_func_grad(iterate),
                                 penalty_parameter,
                                 SLEQP_CAUCHY_OBJECTIVE_TYPE_DEFAULT));

  ASSERT_CALL(sleqp_cauchy_get_working_set(cauchy_data,
                                           iterate));

  ck_assert_int_eq(sleqp_working_set_get_variable_state(working_set, 0), SLEQP_ACTIVE_LOWER);
  ck_assert_int_eq(sleqp_working_set_get_variable_state(working_set, 1), SLEQP_INACTIVE);
}
END_TEST

void working_set_var_teardown()
{
  ASSERT_CALL(sleqp_cauchy_release(&cauchy_data));

  ASSERT_CALL(sleqp_lpi_free(&lp_interface));

  ASSERT_CALL(sleqp_iterate_release(&iterate));

  ASSERT_CALL(sleqp_problem_release(&problem));

  ASSERT_CALL(sleqp_options_release(&options));

  ASSERT_CALL(sleqp_params_release(&params));

  quadfunc_teardown();
}

Suite* working_set_var_test_suite()
{
  Suite *suite;
  TCase *tc_working_set_var;

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
