#include <check.h>
#include <stdlib.h>

#include "aug_jac/standard_aug_jac.h"
#include "cauchy/standard_cauchy.h"
#include "cmp.h"
#include "dual_estimation/dual_estimation_lsq.h"
#include "fact/fact.h"
#include "lp/lpi.h"
#include "mem.h"
#include "util.h"

#include "quadcons_fixture.h"
#include "test_common.h"

SleqpSettings* settings;
SleqpProblem* problem;
SleqpIterate* iterate;
SleqpCauchy* cauchy_data;

SleqpVec* cauchy_direction;

void
constrained_setup()
{
  quadconsfunc_setup();

  ASSERT_CALL(sleqp_settings_create(&settings));

  ASSERT_CALL(sleqp_problem_create_simple(&problem,
                                          quadconsfunc,
                                          quadconsfunc_var_lb,
                                          quadconsfunc_var_ub,
                                          quadconsfunc_cons_lb,
                                          quadconsfunc_cons_ub,
                                          settings));

  ASSERT_CALL(sleqp_iterate_create(&iterate, problem, quadconsfunc_x));

  ASSERT_CALL(
    sleqp_set_and_evaluate(problem, iterate, SLEQP_VALUE_REASON_NONE, NULL));

  ASSERT_CALL(
    sleqp_standard_cauchy_create(&cauchy_data, problem, settings));

  ASSERT_CALL(sleqp_vec_create(&cauchy_direction, 2, 2));
}

START_TEST(test_working_set)
{
  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);

  double trust_radius = 0.1, penalty_parameter = 1.;

  ASSERT_CALL(sleqp_cauchy_set_iterate(cauchy_data, iterate, trust_radius));

  ASSERT_CALL(sleqp_cauchy_solve(cauchy_data,
                                 sleqp_iterate_obj_grad(iterate),
                                 penalty_parameter,
                                 SLEQP_CAUCHY_OBJTYPE_DEFAULT));

  ASSERT_CALL(sleqp_cauchy_working_set(cauchy_data, iterate));

  ck_assert_int_eq(sleqp_working_set_cons_state(working_set, 0),
                   SLEQP_INACTIVE);
  ck_assert_int_eq(sleqp_working_set_cons_state(working_set, 1),
                   SLEQP_ACTIVE_UPPER);
}
END_TEST

START_TEST(test_dual_variable)
{
  SleqpVec* cons_dual = sleqp_iterate_cons_dual(iterate);

  SleqpFact* fact;
  SleqpAugJac* jacobian;
  SleqpDualEstimation* estimation_data;

  double trust_radius = 0.1, penalty_parameter = 1.;

  ASSERT_CALL(sleqp_cauchy_set_iterate(cauchy_data, iterate, trust_radius));

  ASSERT_CALL(sleqp_cauchy_solve(cauchy_data,
                                 sleqp_iterate_obj_grad(iterate),
                                 penalty_parameter,
                                 SLEQP_CAUCHY_OBJTYPE_DEFAULT));

  ASSERT_CALL(sleqp_cauchy_working_set(cauchy_data, iterate));

  ASSERT_CALL(sleqp_cauchy_lp_step(cauchy_data, cauchy_direction));

  ASSERT_CALL(sleqp_fact_create_default(&fact, settings));

  ASSERT_CALL(sleqp_standard_aug_jac_create(&jacobian, problem, settings, fact));

  ASSERT_CALL(sleqp_aug_jac_set_iterate(jacobian, iterate));

  ASSERT_CALL(
    sleqp_dual_estimation_lsq_create(&estimation_data, problem, jacobian));

  ASSERT_CALL(sleqp_estimate_duals(estimation_data,
                                   iterate,
                                   sleqp_iterate_cons_dual(iterate),
                                   sleqp_iterate_vars_dual(iterate)));

  ck_assert_int_eq(cons_dual->dim, 2);
  ck_assert_int_eq(cons_dual->nnz, 1);

  ck_assert_int_eq(cons_dual->dim, 2);

  ck_assert(sleqp_is_eq(cons_dual->data[0], 0.4142135623, 1e-8));

  ASSERT_CALL(sleqp_dual_estimation_release(&estimation_data));

  ASSERT_CALL(sleqp_aug_jac_release(&jacobian));

  ASSERT_CALL(sleqp_fact_release(&fact));
}
END_TEST

void
constrained_teardown()
{
  ASSERT_CALL(sleqp_vec_free(&cauchy_direction));

  ASSERT_CALL(sleqp_cauchy_release(&cauchy_data));

  ASSERT_CALL(sleqp_iterate_release(&iterate));

  ASSERT_CALL(sleqp_problem_release(&problem));

  ASSERT_CALL(sleqp_settings_release(&settings));

  quadconsfunc_teardown();
}

Suite*
constrained_cauchy_test_suite()
{
  Suite* suite;
  TCase* tc_cons;

  suite = suite_create("Constrained tests");

  tc_cons = tcase_create("Constrained Cauchy direction");

  tcase_add_checked_fixture(tc_cons, constrained_setup, constrained_teardown);

  tcase_add_test(tc_cons, test_working_set);
  tcase_add_test(tc_cons, test_dual_variable);

  suite_add_tcase(suite, tc_cons);

  return suite;
}

TEST_MAIN(constrained_cauchy_test_suite)
