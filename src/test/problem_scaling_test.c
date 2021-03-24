#include <stdlib.h>
#include <check.h>
#include <fenv.h>

#include "sleqp_problem_scaling.h"
#include "sleqp_deriv_check.h"
#include "sleqp_util.h"

#include "lp/sleqp_lpi.h"

#include "test_common.h"

#include "quadcons_fixture.h"

SleqpParams* params;
SleqpOptions* options;

SleqpScalingData* scaling;
SleqpProblemScaling* problem_scaling;
SleqpProblem* scaled_problem;

SleqpProblem* problem;
SleqpIterate* iterate;

void problem_scaling_setup()
{
  quadconsfunc_setup();

  ASSERT_CALL(sleqp_params_create(&params));

  ASSERT_CALL(sleqp_options_create(&options));

  ASSERT_CALL(sleqp_problem_create(&problem,
                                   quadconsfunc,
                                   quadconsfunc_var_lb,
                                   quadconsfunc_var_ub,
                                   quadconsfunc_cons_lb,
                                   quadconsfunc_cons_ub));

  ASSERT_CALL(sleqp_scaling_create(&scaling,
                                   problem->num_variables,
                                   problem->num_constraints));

  ASSERT_CALL(sleqp_scaling_set_func_weight(scaling,
                                            2));

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

  ASSERT_CALL(sleqp_iterate_create(&iterate,
                                   problem,
                                   quadconsfunc_x));

  ASSERT_CALL(sleqp_set_and_evaluate(problem,
                                     iterate,
                                     SLEQP_VALUE_REASON_INIT));
}

START_TEST(test_overflow)
{
  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling, 0, 10000));

  SleqpSparseVec* point;

  ASSERT_CALL(sleqp_sparse_vector_create(&point, 2, 2));

  ASSERT_CALL(sleqp_sparse_vector_push(point, 0, 1.));

  SleqpFunc* func = scaled_problem->func;

  int func_grad_nnz, cons_val_nnz, cons_jac_nnz;

  SLEQP_RETCODE retcode = sleqp_func_set_value(func,
                                               point,
                                               SLEQP_VALUE_REASON_NONE,
                                               &func_grad_nnz,
                                               &cons_val_nnz,
                                               &cons_jac_nnz);

  ASSERT_CALL(sleqp_sparse_vector_free(&point));

  ck_assert_int_eq(retcode, SLEQP_MATH_ERROR);
}
END_TEST

START_TEST(test_underflow_warning)
{
  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling, 0, -10000));

  SleqpSparseVec* point;

  ASSERT_CALL(sleqp_sparse_vector_create(&point, 2, 2));

  ASSERT_CALL(sleqp_sparse_vector_push(point, 0, 1.));

  SleqpFunc* func = scaled_problem->func;

  int func_grad_nnz, cons_val_nnz, cons_jac_nnz;

  SLEQP_RETCODE retcode = sleqp_func_set_value(func,
                                               point,
                                               SLEQP_VALUE_REASON_NONE,
                                               &func_grad_nnz,
                                               &cons_val_nnz,
                                               &cons_jac_nnz);

  ASSERT_CALL(sleqp_sparse_vector_free(&point));

  ck_assert_int_eq(retcode, SLEQP_OKAY);
}
END_TEST

START_TEST(test_underflow_error)
{
  ASSERT_CALL(sleqp_options_set_int(options,
                                    SLEQP_OPTION_INT_FLOAT_ERROR_FLAGS,
                                    FE_ALL_EXCEPT));

  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling, 0, -10000));

  SleqpSparseVec* point;

  ASSERT_CALL(sleqp_sparse_vector_create(&point, 2, 2));

  ASSERT_CALL(sleqp_sparse_vector_push(point, 0, 1.));

  SleqpFunc* func = scaled_problem->func;

  int func_grad_nnz, cons_val_nnz, cons_jac_nnz;

  SLEQP_RETCODE retcode = sleqp_func_set_value(func,
                                               point,
                                               SLEQP_VALUE_REASON_NONE,
                                               &func_grad_nnz,
                                               &cons_val_nnz,
                                               &cons_jac_nnz);

  ASSERT_CALL(sleqp_sparse_vector_free(&point));

  ck_assert_int_eq(retcode, SLEQP_MATH_ERROR);
}
END_TEST

START_TEST(test_first_order_deriv)
{
  SleqpIterate* scaled_iterate;

  SleqpDerivCheckData* deriv_check_data;

  ASSERT_CALL(sleqp_iterate_create(&scaled_iterate,
                                   scaled_problem,
                                   quadconsfunc_x));

  ASSERT_CALL(sleqp_scale_point(scaling,
                                sleqp_iterate_get_primal(scaled_iterate)));

  ASSERT_CALL(sleqp_set_and_evaluate(scaled_problem,
                                     scaled_iterate,
                                     SLEQP_VALUE_REASON_NONE));

  ASSERT_CALL(sleqp_deriv_checker_create(&deriv_check_data,
                                         scaled_problem,
                                         params));

  ASSERT_CALL(sleqp_deriv_check_first_order(deriv_check_data,
                                            scaled_iterate));

  ASSERT_CALL(sleqp_deriv_checker_free(&deriv_check_data));

  ASSERT_CALL(sleqp_iterate_release(&scaled_iterate));
}
END_TEST

START_TEST(test_second_order_deriv)
{
  SleqpIterate* scaled_iterate;

  SleqpDerivCheckData* deriv_check_data;

  ASSERT_CALL(sleqp_iterate_create(&scaled_iterate,
                                   scaled_problem,
                                   quadconsfunc_x));

  ASSERT_CALL(sleqp_scale_point(scaling,
                                sleqp_iterate_get_primal(scaled_iterate)));

  ASSERT_CALL(sleqp_set_and_evaluate(scaled_problem,
                                     scaled_iterate,
                                     SLEQP_VALUE_REASON_NONE));

  ASSERT_CALL(sleqp_deriv_checker_create(&deriv_check_data,
                                         scaled_problem,
                                         params));

  ASSERT_CALL(sleqp_deriv_check_second_order_exhaustive(deriv_check_data,
                                                        scaled_iterate));

  ASSERT_CALL(sleqp_deriv_checker_free(&deriv_check_data));

  ASSERT_CALL(sleqp_iterate_release(&scaled_iterate));
}
END_TEST

void problem_scaling_teardown()
{
  ASSERT_CALL(sleqp_iterate_release(&iterate));

  ASSERT_CALL(sleqp_problem_scaling_release(&problem_scaling));

  ASSERT_CALL(sleqp_scaling_release(&scaling));

  ASSERT_CALL(sleqp_problem_free(&problem));

  ASSERT_CALL(sleqp_options_release(&options));

  ASSERT_CALL(sleqp_params_release(&params));

  quadconsfunc_teardown();
}

Suite* problem_scaling_test_suite()
{
  Suite *suite;
  TCase *tc_scale_invalid;
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
