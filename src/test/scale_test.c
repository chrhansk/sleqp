#include <stdlib.h>
#include <check.h>

#include "sleqp.h"

#include "lp/sleqp_lpi.h"

#include "test_common.h"

#include "quadcons_fixture.h"

SleqpParams* params;

SleqpScalingData* scaling;
SleqpProblem* scaled_problem;

SleqpProblem* problem;
SleqpIterate* iterate;

void scaling_setup()
{
  quadconsfunc_setup();

  ASSERT_CALL(sleqp_params_create(&params));

  ASSERT_CALL(sleqp_problem_create(&problem,
                                   quadconsfunc,
                                   params,
                                   quadconsfunc_var_lb,
                                   quadconsfunc_var_ub,
                                   quadconsfunc_cons_lb,
                                   quadconsfunc_cons_ub));

  ASSERT_CALL(sleqp_scaling_create(&scaling,
                                   problem,
                                   params));

  ASSERT_CALL(sleqp_scaling_set_func_weight(scaling,
                                            2));

  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling, 0, -1));
  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling, 1, -6));

  ASSERT_CALL(sleqp_scaling_set_cons_weight(scaling, 0, 7));
  ASSERT_CALL(sleqp_scaling_set_cons_weight(scaling, 1, -1));

  ASSERT_CALL(sleqp_scaling_flush(scaling));

  scaled_problem = sleqp_scaling_get_scaled_problem(scaling);

  ASSERT_CALL(sleqp_iterate_create(&iterate,
                                   problem,
                                   quadconsfunc_x));

  ASSERT_CALL(sleqp_set_and_evaluate(problem,
                                     iterate,
                                     SLEQP_VALUE_REASON_INIT));
}

START_TEST(test_func_grad_invalid)
{
  SleqpSparseVec* func_grad;

  ASSERT_CALL(sleqp_sparse_vector_create(&func_grad, 2, 2));

  ASSERT_CALL(sleqp_sparse_vector_copy(iterate->func_grad, func_grad));

  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling, 0, 10000));

  SLEQP_RETCODE scale_retcode = sleqp_scale_func_grad(scaling, func_grad);

  ck_assert_int_eq(scale_retcode, SLEQP_MATH_ERROR);

  ASSERT_CALL(sleqp_sparse_vector_free(&func_grad));
}
END_TEST

START_TEST(test_func_val_inverse)
{
  double func_val = iterate->func_val;

  double scaled_func_val = sleqp_scale_func_val(scaling, func_val);

  double unscaled_func_val = sleqp_unscale_func_val(scaling, scaled_func_val);

  ck_assert(func_val == unscaled_func_val);

}
END_TEST

START_TEST(test_func_grad_inverse)
{
  SleqpSparseVec* func_grad;

  ASSERT_CALL(sleqp_sparse_vector_create(&func_grad, 2, 2));

  ASSERT_CALL(sleqp_sparse_vector_copy(iterate->func_grad, func_grad));


  ASSERT_CALL(sleqp_scale_func_grad(scaling,
                                    func_grad));

  ASSERT_CALL(sleqp_unscale_func_grad(scaling,
                                      func_grad));

  ck_assert(sleqp_sparse_vector_eq(iterate->func_grad,
                                   func_grad,
                                   0.));

  ASSERT_CALL(sleqp_sparse_vector_free(&func_grad));
}
END_TEST

START_TEST(test_cons_val_inverse)
{
  SleqpSparseVec* cons_val;

  ASSERT_CALL(sleqp_sparse_vector_create(&cons_val, 2, 2));

  ASSERT_CALL(sleqp_sparse_vector_copy(iterate->cons_val, cons_val));


  ASSERT_CALL(sleqp_scale_cons_val(scaling,
                                    cons_val));

  ASSERT_CALL(sleqp_unscale_cons_val(scaling,
                                      cons_val));

  ck_assert(sleqp_sparse_vector_eq(iterate->cons_val,
                                   cons_val,
                                   0.));

  ASSERT_CALL(sleqp_sparse_vector_free(&cons_val));
}
END_TEST

START_TEST(test_cons_jac_inverse)
{
  SleqpSparseMatrix* cons_jac;

  ASSERT_CALL(sleqp_sparse_matrix_create(&cons_jac,
                                         2,
                                         2,
                                         4));

  ASSERT_CALL(sleqp_sparse_matrix_copy(iterate->cons_jac,
                                       cons_jac));

  ASSERT_CALL(sleqp_scale_cons_jac(scaling, cons_jac));

  ASSERT_CALL(sleqp_unscale_cons_jac(scaling, cons_jac));

  ck_assert(sleqp_sparse_matrix_eq(iterate->cons_jac,
                                   cons_jac,
                                   0.));

  ASSERT_CALL(sleqp_sparse_matrix_free(&cons_jac));


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
                                scaled_iterate->primal));

  ASSERT_CALL(sleqp_set_and_evaluate(scaled_problem,
                                     scaled_iterate,
                                     SLEQP_VALUE_REASON_NONE));

  ASSERT_CALL(sleqp_deriv_checker_create(&deriv_check_data,
                                         scaled_problem,
                                         params));

  ASSERT_CALL(sleqp_deriv_check_first_order(deriv_check_data,
                                            scaled_iterate));

  ASSERT_CALL(sleqp_deriv_checker_free(&deriv_check_data));

  ASSERT_CALL(sleqp_iterate_free(&scaled_iterate));
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
                                scaled_iterate->primal));

  ASSERT_CALL(sleqp_set_and_evaluate(scaled_problem,
                                     scaled_iterate,
                                     SLEQP_VALUE_REASON_NONE));

  ASSERT_CALL(sleqp_deriv_checker_create(&deriv_check_data,
                                         scaled_problem,
                                         params));

  ASSERT_CALL(sleqp_deriv_check_second_order(deriv_check_data,
                                             scaled_iterate));

  ASSERT_CALL(sleqp_deriv_checker_free(&deriv_check_data));

  ASSERT_CALL(sleqp_iterate_free(&scaled_iterate));
}
END_TEST

void scaling_teardown()
{
  ASSERT_CALL(sleqp_iterate_free(&iterate));

  ASSERT_CALL(sleqp_scaling_free(&scaling));

  ASSERT_CALL(sleqp_problem_free(&problem));

  ASSERT_CALL(sleqp_params_free(&params));

  quadconsfunc_teardown();
}

Suite* scaling_test_suite()
{
  Suite *suite;
  TCase *tc_scale_inv;
  TCase* tc_scale_deriv;

  suite = suite_create("Scaling tests");

  tc_scale_inv = tcase_create("Scaling inverse tests");

  tc_scale_deriv = tcase_create("Scaling derivative tests");

  tcase_add_checked_fixture(tc_scale_inv,
                            scaling_setup,
                            scaling_teardown);

  tcase_add_checked_fixture(tc_scale_deriv,
                            scaling_setup,
                            scaling_teardown);
  tcase_add_test(tc_scale_inv, test_func_val_inverse);
  tcase_add_test(tc_scale_inv, test_func_grad_invalid);
  tcase_add_test(tc_scale_inv, test_func_grad_inverse);
  tcase_add_test(tc_scale_inv, test_cons_val_inverse);
  tcase_add_test(tc_scale_inv, test_cons_jac_inverse);

  tcase_add_test(tc_scale_deriv, test_first_order_deriv);
  tcase_add_test(tc_scale_deriv, test_second_order_deriv);

  suite_add_tcase(suite, tc_scale_inv);

  suite_add_tcase(suite, tc_scale_deriv);

  return suite;
}

int main()
{
  int num_fails;
  Suite* suite;
  SRunner* srunner;

  suite = scaling_test_suite();
  srunner = srunner_create(suite);

  srunner_set_fork_status(srunner, CK_NOFORK);
  srunner_run_all(srunner, CK_NORMAL);

  num_fails = srunner_ntests_failed(srunner);

  srunner_free(srunner);

  return (num_fails > 0) ? EXIT_FAILURE : EXIT_SUCCESS;
}
