#include <check.h>
#include <stdlib.h>

#include "aug_jac/standard_aug_jac.h"
#include "cauchy/standard_cauchy.h"
#include "cmp.h"
#include "dual_estimation/dual_estimation_lsq.h"
#include "factorization/factorization.h"
#include "mem.h"
#include "test_common.h"
#include "util.h"

#include "quadfunc_fixture.h"

START_TEST(test_simply_constrained_dual_estimation)
{
  SleqpParams* params;
  SleqpOptions* options;
  SleqpProblem* problem;
  SleqpIterate* iterate;
  SleqpCauchy* cauchy_data;
  SleqpWorkingSet* working_set;
  SleqpFact* factorization;
  SleqpAugJac* aug_jac;

  SleqpDualEstimation* estimation_data;

  double penalty_parameter = 1., trust_radius = 0.1;

  ASSERT_CALL(sleqp_params_create(&params));

  ASSERT_CALL(sleqp_options_create(&options));

  ASSERT_CALL(sleqp_problem_create_simple(&problem,
                                          quadfunc,
                                          params,
                                          quadfunc_var_lb,
                                          quadfunc_var_ub,
                                          quadfunc_cons_lb,
                                          quadfunc_cons_ub));

  ASSERT_CALL(sleqp_iterate_create(&iterate, problem, quadfunc_x));

  ASSERT_CALL(
    sleqp_set_and_evaluate(problem, iterate, SLEQP_VALUE_REASON_NONE, NULL));

  ASSERT_CALL(
    sleqp_standard_cauchy_create(&cauchy_data, problem, params, options));

  ASSERT_CALL(sleqp_working_set_create(&working_set, problem));

  ASSERT_CALL(sleqp_fact_create_default(&factorization, params));

  ASSERT_CALL(
    sleqp_standard_aug_jac_create(&aug_jac, problem, params, factorization));

  ASSERT_CALL(
    sleqp_dual_estimation_lsq_create(&estimation_data, problem, aug_jac));

  ASSERT_CALL(sleqp_cauchy_set_iterate(cauchy_data, iterate, trust_radius));

  ASSERT_CALL(sleqp_cauchy_solve(cauchy_data,
                                 sleqp_iterate_obj_grad(iterate),
                                 penalty_parameter,
                                 SLEQP_CAUCHY_OBJECTIVE_TYPE_DEFAULT));

  ASSERT_CALL(sleqp_cauchy_working_set(cauchy_data, iterate));

  ASSERT_CALL(sleqp_aug_jac_set_iterate(aug_jac, iterate));

  /*
    ASSERT_CALL(
      sleqp_dual_estimation_compute(estimation_data, iterate, NULL, aug_jac));
    */

  ASSERT_CALL(sleqp_estimate_duals(estimation_data,
                                   iterate,
                                   sleqp_iterate_cons_dual(iterate),
                                   sleqp_iterate_vars_dual(iterate)));

  SleqpVec* vars_dual = sleqp_iterate_vars_dual(iterate);

  ck_assert(sleqp_vec_at(vars_dual, 0));
  ck_assert(sleqp_vec_at(vars_dual, 1));

  double tolerance = 1e-8;

  ck_assert(sleqp_is_eq(*sleqp_vec_at(vars_dual, 0), -2., tolerance));
  ck_assert(sleqp_is_eq(*sleqp_vec_at(vars_dual, 1), -4., tolerance));

  ASSERT_CALL(sleqp_dual_estimation_release(&estimation_data));

  ASSERT_CALL(sleqp_aug_jac_release(&aug_jac));

  ASSERT_CALL(sleqp_fact_release(&factorization));

  ASSERT_CALL(sleqp_working_set_release(&working_set));

  ASSERT_CALL(sleqp_cauchy_release(&cauchy_data));

  ASSERT_CALL(sleqp_iterate_release(&iterate));

  ASSERT_CALL(sleqp_problem_release(&problem));

  ASSERT_CALL(sleqp_options_release(&options));

  ASSERT_CALL(sleqp_params_release(&params));
}
END_TEST

Suite*
dual_estimation_test_suite()
{
  Suite* suite;
  TCase* tc_dual_estimation;

  suite = suite_create("Dual estimation tests");

  tc_dual_estimation = tcase_create("Simply constrained");

  tcase_add_checked_fixture(tc_dual_estimation,
                            quadfunc_setup,
                            quadfunc_teardown);

  tcase_add_test(tc_dual_estimation, test_simply_constrained_dual_estimation);
  suite_add_tcase(suite, tc_dual_estimation);

  return suite;
}

TEST_MAIN(dual_estimation_test_suite)
