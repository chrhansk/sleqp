#include <check.h>
#include <stdlib.h>

#include "lp/lpi.h"
#include "scale.h"
#include "util.h"

#include "test_common.h"

#include "quadcons_fixture.h"

SleqpParams* params;

SleqpScaling* scaling;

SleqpProblem* problem;
SleqpIterate* iterate;

void
scaling_setup()
{
  quadconsfunc_setup();

  ASSERT_CALL(sleqp_params_create(&params));

  ASSERT_CALL(sleqp_problem_create_simple(&problem,
                                          quadconsfunc,
                                          params,
                                          quadconsfunc_var_lb,
                                          quadconsfunc_var_ub,
                                          quadconsfunc_cons_lb,
                                          quadconsfunc_cons_ub));

  const int num_variables   = sleqp_problem_num_vars(problem);
  const int num_constraints = sleqp_problem_num_cons(problem);

  ASSERT_CALL(sleqp_scaling_create(&scaling, num_variables, num_constraints));

  ASSERT_CALL(sleqp_scaling_set_obj_weight(scaling, 2));

  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling, 0, -1));
  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling, 1, -6));

  ASSERT_CALL(sleqp_scaling_set_cons_weight(scaling, 0, 7));
  ASSERT_CALL(sleqp_scaling_set_cons_weight(scaling, 1, -1));

  ASSERT_CALL(sleqp_iterate_create(&iterate, problem, quadconsfunc_x));

  ASSERT_CALL(
    sleqp_set_and_evaluate(problem, iterate, SLEQP_VALUE_REASON_INIT, NULL));
}

START_TEST(test_nominal_large)
{
  double nominal_values[] = {7.9, 15.9};

  ASSERT_CALL(
    sleqp_scaling_set_var_weights_from_nominal(scaling, nominal_values));

  int* var_weights = sleqp_scaling_var_weights(scaling);

  ck_assert_int_eq(var_weights[0], 3);
  ck_assert_int_eq(var_weights[1], 4);
}
END_TEST

START_TEST(test_nominal_scale)
{
  const double eps        = sleqp_params_value(params, SLEQP_PARAM_EPS);
  double nominal_values[] = {4.1, 15.9};

  ASSERT_CALL(
    sleqp_scaling_set_var_weights_from_nominal(scaling, nominal_values));

  SleqpVec* primal = sleqp_iterate_primal(iterate);

  ASSERT_CALL(sleqp_vec_clear(primal));

  ASSERT_CALL(sleqp_vec_push(primal, 0, nominal_values[0] - 1e-5));

  ASSERT_CALL(sleqp_vec_push(primal, 1, nominal_values[1] - 1e-5));

  ASSERT_CALL(sleqp_scale_iterate(scaling, iterate, false));

  for (int index = 0; index < primal->nnz; ++index)
  {
    ck_assert(sleqp_is_geq(primal->data[index], .5, eps));
    ck_assert(sleqp_is_leq(primal->data[index], 1., eps));
  }
}
END_TEST

START_TEST(test_nominal_neg)
{
  double nominal_values[] = {-7.9, -15.9};

  ASSERT_CALL(
    sleqp_scaling_set_var_weights_from_nominal(scaling, nominal_values));

  int* var_weights = sleqp_scaling_var_weights(scaling);

  ck_assert_int_eq(var_weights[0], 3);
  ck_assert_int_eq(var_weights[1], 4);
}
END_TEST

START_TEST(test_nominal_small)
{
  double nominal_values[] = {0.51, 0.26};

  ASSERT_CALL(
    sleqp_scaling_set_var_weights_from_nominal(scaling, nominal_values));

  int* var_weights = sleqp_scaling_var_weights(scaling);

  ck_assert_int_eq(var_weights[0], 0);
  ck_assert_int_eq(var_weights[1], -1);
}
END_TEST

START_TEST(test_obj_val_inverse)
{
  double obj_val = sleqp_iterate_obj_val(iterate);

  double scaled_obj_val = sleqp_scale_obj_val(scaling, obj_val);

  double unscaled_obj_val = sleqp_unscale_obj_val(scaling, scaled_obj_val);

  ck_assert(obj_val == unscaled_obj_val);
}
END_TEST

START_TEST(test_nominal_scale_obj_val)
{
  double nominal_obj_val = 17.;

  const double eps = sleqp_params_value(params, SLEQP_PARAM_EPS);

  ASSERT_CALL(
    sleqp_scaling_set_obj_weight_from_nominal(scaling, nominal_obj_val));

  double scaled_obj_val = sleqp_scale_obj_val(scaling, nominal_obj_val);

  ck_assert(sleqp_is_leq(scaled_obj_val, 1., eps));
  ck_assert(sleqp_is_geq(scaled_obj_val, .5, eps));
}
END_TEST

START_TEST(test_nominal_scale_cons_vals)
{
  double nominal_cons_vals[] = {17., .4};

  const int num_constraints = 2.;
  const double eps          = sleqp_params_value(params, SLEQP_PARAM_EPS);

  SleqpVec* cons_vals;

  ASSERT_CALL(sleqp_vec_create(&cons_vals, num_constraints, num_constraints));

  ASSERT_CALL(
    sleqp_vec_from_raw(cons_vals, nominal_cons_vals, num_constraints, eps));

  ASSERT_CALL(
    sleqp_scaling_set_cons_weights_from_nominal(scaling, nominal_cons_vals));

  ASSERT_CALL(sleqp_scale_cons_val(scaling, cons_vals));

  for (int k = 0; k < cons_vals->nnz; ++k)
  {
    const double value = cons_vals->data[k];
    ck_assert(sleqp_is_leq(value, 1., eps));
    ck_assert(sleqp_is_geq(value, .5, eps));
  }

  ASSERT_CALL(sleqp_vec_free(&cons_vals));
}
END_TEST

START_TEST(test_obj_grad_inverse)
{
  SleqpVec* obj_grad;

  ASSERT_CALL(sleqp_vec_create(&obj_grad, 2, 2));

  ASSERT_CALL(sleqp_vec_copy(sleqp_iterate_obj_grad(iterate), obj_grad));

  ASSERT_CALL(sleqp_scale_obj_grad(scaling, obj_grad));

  ASSERT_CALL(sleqp_unscale_obj_grad(scaling, obj_grad));

  ck_assert(sleqp_vec_eq(sleqp_iterate_obj_grad(iterate), obj_grad, 0.));

  ASSERT_CALL(sleqp_vec_free(&obj_grad));
}
END_TEST

START_TEST(test_cons_val_inverse)
{
  SleqpVec* cons_val;

  ASSERT_CALL(sleqp_vec_create(&cons_val, 2, 2));

  ASSERT_CALL(sleqp_vec_copy(sleqp_iterate_cons_val(iterate), cons_val));

  ASSERT_CALL(sleqp_scale_cons_val(scaling, cons_val));

  ASSERT_CALL(sleqp_unscale_cons_val(scaling, cons_val));

  ck_assert(sleqp_vec_eq(sleqp_iterate_cons_val(iterate), cons_val, 0.));

  ASSERT_CALL(sleqp_vec_free(&cons_val));
}
END_TEST

START_TEST(test_cons_jac_inverse)
{
  SleqpSparseMatrix* cons_jac;

  ASSERT_CALL(sleqp_sparse_matrix_create(&cons_jac, 2, 2, 4));

  ASSERT_CALL(
    sleqp_sparse_matrix_copy(sleqp_iterate_cons_jac(iterate), cons_jac));

  ASSERT_CALL(sleqp_scale_cons_jac(scaling, cons_jac));

  ASSERT_CALL(sleqp_unscale_cons_jac(scaling, cons_jac));

  ck_assert(
    sleqp_sparse_matrix_eq(sleqp_iterate_cons_jac(iterate), cons_jac, 0.));

  ASSERT_CALL(sleqp_sparse_matrix_release(&cons_jac));
}
END_TEST

void
scaling_teardown()
{
  ASSERT_CALL(sleqp_iterate_release(&iterate));

  ASSERT_CALL(sleqp_scaling_release(&scaling));

  ASSERT_CALL(sleqp_problem_release(&problem));

  ASSERT_CALL(sleqp_params_release(&params));

  quadconsfunc_teardown();
}

Suite*
scaling_test_suite()
{
  Suite* suite;
  TCase* tc_nominal;
  TCase* tc_scale_inv;

  suite = suite_create("Scaling tests");

  tc_nominal = tcase_create("Nominal values tests");

  tc_scale_inv = tcase_create("Scaling inverse tests");

  tcase_add_checked_fixture(tc_nominal, scaling_setup, scaling_teardown);

  tcase_add_checked_fixture(tc_scale_inv, scaling_setup, scaling_teardown);

  tcase_add_test(tc_nominal, test_nominal_large);
  tcase_add_test(tc_nominal, test_nominal_scale);
  tcase_add_test(tc_nominal, test_nominal_neg);
  tcase_add_test(tc_nominal, test_nominal_small);
  tcase_add_test(tc_nominal, test_nominal_scale_obj_val);
  tcase_add_test(tc_nominal, test_nominal_scale_cons_vals);

  tcase_add_test(tc_scale_inv, test_obj_val_inverse);
  tcase_add_test(tc_scale_inv, test_obj_grad_inverse);
  tcase_add_test(tc_scale_inv, test_cons_val_inverse);
  tcase_add_test(tc_scale_inv, test_cons_jac_inverse);

  suite_add_tcase(suite, tc_nominal);

  suite_add_tcase(suite, tc_scale_inv);

  return suite;
}

TEST_MAIN(scaling_test_suite)
