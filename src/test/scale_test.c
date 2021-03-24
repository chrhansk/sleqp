#include <stdlib.h>
#include <check.h>

#include "sleqp_util.h"

#include "lp/sleqp_lpi.h"

#include "test_common.h"

#include "quadcons_fixture.h"

SleqpParams* params;

SleqpScalingData* scaling;

SleqpProblem* problem;
SleqpIterate* iterate;

void scaling_setup()
{
  quadconsfunc_setup();

  ASSERT_CALL(sleqp_params_create(&params));

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

  ASSERT_CALL(sleqp_iterate_create(&iterate,
                                   problem,
                                   quadconsfunc_x));

  ASSERT_CALL(sleqp_set_and_evaluate(problem,
                                     iterate,
                                     SLEQP_VALUE_REASON_INIT));
}

START_TEST(test_nominal_large)
{
  double nominal_values[] = {7.9, 15.9};

  ASSERT_CALL(sleqp_scaling_set_var_weights_from_nominal(scaling,
                                                         nominal_values));

  int* var_weights = sleqp_scaling_get_var_weights(scaling);

  ck_assert_int_eq(var_weights[0], 3);
  ck_assert_int_eq(var_weights[1], 4);
}
END_TEST

START_TEST(test_nominal_scale)
{
  const double eps = sleqp_params_get(params, SLEQP_PARAM_EPS);
  double nominal_values[] = {4.1, 15.9};

  ASSERT_CALL(sleqp_scaling_set_var_weights_from_nominal(scaling,
                                                         nominal_values));

  SleqpSparseVec* primal = sleqp_iterate_get_primal(iterate);

  ASSERT_CALL(sleqp_sparse_vector_clear(primal));

  ASSERT_CALL(sleqp_sparse_vector_push(primal,
                                       0,
                                       nominal_values[0] - 1e-5));

  ASSERT_CALL(sleqp_sparse_vector_push(primal,
                                       1,
                                       nominal_values[1] - 1e-5));

  ASSERT_CALL(sleqp_scale_iterate(scaling, iterate));

  for(int index = 0; index < primal->nnz; ++index)
  {
    ck_assert(sleqp_is_geq(primal->data[index], .5, eps));
    ck_assert(sleqp_is_leq(primal->data[index], 1., eps));
  }


}
END_TEST

START_TEST(test_nominal_neg)
{
  double nominal_values[] = {-7.9, -15.9};

  ASSERT_CALL(sleqp_scaling_set_var_weights_from_nominal(scaling,
                                                         nominal_values));

  int* var_weights = sleqp_scaling_get_var_weights(scaling);

  ck_assert_int_eq(var_weights[0], 3);
  ck_assert_int_eq(var_weights[1], 4);
}
END_TEST

START_TEST(test_nominal_small)
{
  double nominal_values[] = {0.51, 0.26};

  ASSERT_CALL(sleqp_scaling_set_var_weights_from_nominal(scaling,
                                                         nominal_values));

  int* var_weights = sleqp_scaling_get_var_weights(scaling);

  ck_assert_int_eq(var_weights[0], 0);
  ck_assert_int_eq(var_weights[1], -1);
}
END_TEST

START_TEST(test_func_val_inverse)
{
  double func_val = sleqp_iterate_get_func_val(iterate);

  double scaled_func_val = sleqp_scale_func_val(scaling, func_val);

  double unscaled_func_val = sleqp_unscale_func_val(scaling, scaled_func_val);

  ck_assert(func_val == unscaled_func_val);

}
END_TEST

START_TEST(test_nominal_scale_func_val)
{
  double nominal_func_val = 17.;

  const double eps = sleqp_params_get(params,
                                      SLEQP_PARAM_EPS);

  ASSERT_CALL(sleqp_scaling_set_func_weight_from_nominal(scaling,
                                                         nominal_func_val));

  double scaled_func_val = sleqp_scale_func_val(scaling,
                                                nominal_func_val);

  ck_assert(sleqp_is_leq(scaled_func_val, 1., eps));
  ck_assert(sleqp_is_geq(scaled_func_val, .5, eps));
}
END_TEST

START_TEST(test_nominal_scale_cons_vals)
{
  double nominal_cons_vals[] = {17., .4};

  const int num_constraints = 2.;
  const double eps = sleqp_params_get(params,
                                      SLEQP_PARAM_EPS);

  SleqpSparseVec* cons_vals;

  ASSERT_CALL(sleqp_sparse_vector_create(&cons_vals,
                                         num_constraints,
                                         num_constraints));

  ASSERT_CALL(sleqp_sparse_vector_from_raw(cons_vals,
                                           nominal_cons_vals,
                                           num_constraints,
                                           eps));

  ASSERT_CALL(sleqp_scaling_set_cons_weights_from_nominal(scaling,
                                                          nominal_cons_vals));

  ASSERT_CALL(sleqp_scale_cons_val(scaling, cons_vals));

  for(int k = 0; k < cons_vals->nnz; ++k)
  {
    const double value = cons_vals->data[k];
    ck_assert(sleqp_is_leq(value, 1., eps));
    ck_assert(sleqp_is_geq(value, .5, eps));
  }

  ASSERT_CALL(sleqp_sparse_vector_free(&cons_vals));
}
END_TEST

START_TEST(test_func_grad_inverse)
{
  SleqpSparseVec* func_grad;

  ASSERT_CALL(sleqp_sparse_vector_create(&func_grad, 2, 2));

  ASSERT_CALL(sleqp_sparse_vector_copy(sleqp_iterate_get_func_grad(iterate),
                                       func_grad));


  ASSERT_CALL(sleqp_scale_func_grad(scaling,
                                    func_grad));

  ASSERT_CALL(sleqp_unscale_func_grad(scaling,
                                      func_grad));

  ck_assert(sleqp_sparse_vector_eq(sleqp_iterate_get_func_grad(iterate),
                                   func_grad,
                                   0.));

  ASSERT_CALL(sleqp_sparse_vector_free(&func_grad));
}
END_TEST

START_TEST(test_cons_val_inverse)
{
  SleqpSparseVec* cons_val;

  ASSERT_CALL(sleqp_sparse_vector_create(&cons_val, 2, 2));

  ASSERT_CALL(sleqp_sparse_vector_copy(sleqp_iterate_get_cons_val(iterate), cons_val));


  ASSERT_CALL(sleqp_scale_cons_val(scaling,
                                   cons_val));

  ASSERT_CALL(sleqp_unscale_cons_val(scaling,
                                     cons_val));

  ck_assert(sleqp_sparse_vector_eq(sleqp_iterate_get_cons_val(iterate),
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

  ASSERT_CALL(sleqp_sparse_matrix_copy(sleqp_iterate_get_cons_jac(iterate),
                                       cons_jac));

  ASSERT_CALL(sleqp_scale_cons_jac(scaling, cons_jac));

  ASSERT_CALL(sleqp_unscale_cons_jac(scaling, cons_jac));

  ck_assert(sleqp_sparse_matrix_eq(sleqp_iterate_get_cons_jac(iterate),
                                   cons_jac,
                                   0.));

  ASSERT_CALL(sleqp_sparse_matrix_release(&cons_jac));


}
END_TEST

void scaling_teardown()
{
  ASSERT_CALL(sleqp_iterate_release(&iterate));

  ASSERT_CALL(sleqp_scaling_release(&scaling));

  ASSERT_CALL(sleqp_problem_free(&problem));

  ASSERT_CALL(sleqp_params_release(&params));

  quadconsfunc_teardown();
}

Suite* scaling_test_suite()
{
  Suite *suite;
  TCase *tc_nominal;
  TCase *tc_scale_inv;

  suite = suite_create("Scaling tests");

  tc_nominal = tcase_create("Nominal values tests");

  tc_scale_inv = tcase_create("Scaling inverse tests");

  tcase_add_checked_fixture(tc_nominal,
                            scaling_setup,
                            scaling_teardown);

  tcase_add_checked_fixture(tc_scale_inv,
                            scaling_setup,
                            scaling_teardown);

  tcase_add_test(tc_nominal, test_nominal_large);
  tcase_add_test(tc_nominal, test_nominal_scale);
  tcase_add_test(tc_nominal, test_nominal_neg);
  tcase_add_test(tc_nominal, test_nominal_small);
  tcase_add_test(tc_nominal, test_nominal_scale_func_val);
  tcase_add_test(tc_nominal, test_nominal_scale_cons_vals);

  tcase_add_test(tc_scale_inv, test_func_val_inverse);
  tcase_add_test(tc_scale_inv, test_func_grad_inverse);
  tcase_add_test(tc_scale_inv, test_cons_val_inverse);
  tcase_add_test(tc_scale_inv, test_cons_jac_inverse);

  suite_add_tcase(suite, tc_nominal);

  suite_add_tcase(suite, tc_scale_inv);

  return suite;
}

TEST_MAIN(scaling_test_suite)
