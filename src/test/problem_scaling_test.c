#include <check.h>
#include <fenv.h>
#include <stdlib.h>

#include "deriv_check.h"
#include "problem_scaling.h"
#include "util.h"

#include "lp/lpi.h"

#include "test_common.h"

#include "quadcons_fixture.h"

SleqpParams* params;
SleqpOptions* options;

SleqpScaling* scaling;
SleqpProblemScaling* problem_scaling;
SleqpProblem* scaled_problem;

SleqpProblem* problem;
SleqpIterate* iterate;

void
problem_scaling_setup()
{
  quadconsfunc_setup();

  ASSERT_CALL(sleqp_params_create(&params));

  ASSERT_CALL(sleqp_options_create(&options));

  ASSERT_CALL(sleqp_problem_create_simple(&problem,
                                          quadconsfunc,
                                          params,
                                          quadconsfunc_var_lb,
                                          quadconsfunc_var_ub,
                                          quadconsfunc_cons_lb,
                                          quadconsfunc_cons_ub));

  const int num_constraints = sleqp_problem_num_cons(problem);

  ASSERT_CALL(sleqp_scaling_create(&scaling,
                                   sleqp_problem_num_vars(problem),
                                   num_constraints));

  ASSERT_CALL(sleqp_scaling_set_obj_weight(scaling, 2));

  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling, 0, -1));
  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling, 1, -6));

  ASSERT_CALL(sleqp_scaling_set_cons_weight(scaling, 0, 7));
  ASSERT_CALL(sleqp_scaling_set_cons_weight(scaling, 1, -1));

  ASSERT_CALL(sleqp_problem_scaling_create(&problem_scaling,
                                           scaling,
                                           problem,
                                           params,
                                           options));

  ASSERT_CALL(sleqp_problem_scaling_flush(problem_scaling));

  scaled_problem = sleqp_problem_scaling_get_problem(problem_scaling);

  ASSERT_CALL(sleqp_iterate_create(&iterate, problem, quadconsfunc_x));

  ASSERT_CALL(
    sleqp_set_and_evaluate(problem, iterate, SLEQP_VALUE_REASON_INIT, NULL));
}

START_TEST(test_overflow)
{
  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling, 0, 10000));

  SleqpVec* point;

  ASSERT_CALL(sleqp_vec_create(&point, 2, 2));

  ASSERT_CALL(sleqp_vec_push(point, 0, 1.));

  SleqpFunc* func = sleqp_problem_func(scaled_problem);

  bool reject;

  SLEQP_RETCODE retcode
    = sleqp_func_set_value(func, point, SLEQP_VALUE_REASON_NONE, &reject);

  ASSERT_CALL(sleqp_vec_free(&point));

  ck_assert_int_eq(retcode, SLEQP_ERROR);
  ck_assert_int_eq(sleqp_error_type(), SLEQP_MATH_ERROR);
}
END_TEST

START_TEST(test_underflow_warning)
{
  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling, 0, -10000));

  SleqpVec* point;

  ASSERT_CALL(sleqp_vec_create(&point, 2, 2));

  ASSERT_CALL(sleqp_vec_push(point, 0, 1.));

  SleqpFunc* func = sleqp_problem_func(scaled_problem);

  bool reject;

  SLEQP_RETCODE retcode
    = sleqp_func_set_value(func, point, SLEQP_VALUE_REASON_NONE, &reject);

  ASSERT_CALL(sleqp_vec_free(&point));

  ck_assert_int_eq(retcode, SLEQP_OKAY);
}
END_TEST

START_TEST(test_underflow_error)
{
  ASSERT_CALL(sleqp_options_set_enum_value(options,
                                           SLEQP_OPTION_ENUM_FLOAT_ERROR_FLAGS,
                                           FE_ALL_EXCEPT));

  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling, 0, -10000));

  SleqpVec* point;

  ASSERT_CALL(sleqp_vec_create(&point, 2, 2));

  ASSERT_CALL(sleqp_vec_push(point, 0, 1.));

  SleqpFunc* func = sleqp_problem_func(scaled_problem);

  bool reject;

  SLEQP_RETCODE retcode
    = sleqp_func_set_value(func, point, SLEQP_VALUE_REASON_NONE, &reject);

  ASSERT_CALL(sleqp_vec_free(&point));

  ck_assert_int_eq(retcode, SLEQP_ERROR);
  ck_assert_int_eq(sleqp_error_type(), SLEQP_MATH_ERROR);
}
END_TEST

START_TEST(test_first_order_deriv)
{
  SleqpIterate* scaled_iterate;

  SleqpDerivChecker* deriv_check_data;

  ASSERT_CALL(
    sleqp_iterate_create(&scaled_iterate, scaled_problem, quadconsfunc_x));

  ASSERT_CALL(sleqp_scale_point(scaling, sleqp_iterate_primal(scaled_iterate)));

  ASSERT_CALL(sleqp_set_and_evaluate(scaled_problem,
                                     scaled_iterate,
                                     SLEQP_VALUE_REASON_NONE,
                                     NULL));

  ASSERT_CALL(
    sleqp_deriv_checker_create(&deriv_check_data, scaled_problem, params));

  ASSERT_CALL(sleqp_deriv_check_perform(deriv_check_data,
                                        scaled_iterate,
                                        SLEQP_DERIV_CHECK_FIRST));

  ASSERT_CALL(sleqp_deriv_checker_free(&deriv_check_data));

  ASSERT_CALL(sleqp_iterate_release(&scaled_iterate));
}
END_TEST

START_TEST(test_second_order_deriv)
{
  SleqpIterate* scaled_iterate;

  SleqpDerivChecker* deriv_check_data;

  ASSERT_CALL(
    sleqp_iterate_create(&scaled_iterate, scaled_problem, quadconsfunc_x));

  ASSERT_CALL(sleqp_scale_point(scaling, sleqp_iterate_primal(scaled_iterate)));

  ASSERT_CALL(sleqp_set_and_evaluate(scaled_problem,
                                     scaled_iterate,
                                     SLEQP_VALUE_REASON_NONE,
                                     NULL));

  ASSERT_CALL(
    sleqp_deriv_checker_create(&deriv_check_data, scaled_problem, params));

  ASSERT_CALL(sleqp_deriv_check_perform(deriv_check_data,
                                        scaled_iterate,
                                        SLEQP_DERIV_CHECK_SECOND_EXHAUSTIVE));

  ASSERT_CALL(sleqp_deriv_checker_free(&deriv_check_data));

  ASSERT_CALL(sleqp_iterate_release(&scaled_iterate));
}
END_TEST

void
problem_scaling_teardown()
{
  ASSERT_CALL(sleqp_iterate_release(&iterate));

  ASSERT_CALL(sleqp_problem_scaling_release(&problem_scaling));

  ASSERT_CALL(sleqp_scaling_release(&scaling));

  ASSERT_CALL(sleqp_problem_release(&problem));

  ASSERT_CALL(sleqp_options_release(&options));

  ASSERT_CALL(sleqp_params_release(&params));

  quadconsfunc_teardown();
}

Suite*
problem_scaling_test_suite()
{
  Suite* suite;
  TCase* tc_scale_invalid;
  TCase* tc_scale_deriv;

  suite = suite_create("Problem scaling tests");

  tc_scale_deriv = tcase_create("Scaling derivative tests");

  tcase_add_checked_fixture(tc_scale_deriv,
                            problem_scaling_setup,
                            problem_scaling_teardown);

  tc_scale_invalid = tcase_create("Invalid scaling values");

  tcase_add_test(tc_scale_invalid, test_overflow);
  tcase_add_test(tc_scale_invalid, test_underflow_warning);
  tcase_add_test(tc_scale_invalid, test_underflow_error);

  tcase_add_checked_fixture(tc_scale_invalid,
                            problem_scaling_setup,
                            problem_scaling_teardown);

  tcase_add_test(tc_scale_deriv, test_first_order_deriv);
  tcase_add_test(tc_scale_deriv, test_second_order_deriv);

  suite_add_tcase(suite, tc_scale_invalid);
  suite_add_tcase(suite, tc_scale_deriv);

  return suite;
}

TEST_MAIN(problem_scaling_test_suite)
