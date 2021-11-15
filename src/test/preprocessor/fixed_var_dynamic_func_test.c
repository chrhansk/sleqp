#include <check.h>

#include "dyn_rosenbrock_fixture.h"
#include "test_common.h"

#include "cmp.h"
#include "dyn.h"
#include "preprocessor/fixed_var_func.h"

const double accuracy = 1e-8;

int num_fixed         = 1;
int fixed_indices[]   = {0};
double fixed_values[] = {1.};

SleqpFunc* fixed_var_func;

SleqpSparseVec* fixed_initial;

void
setup()
{
  dyn_rosenbrock_setup();

  ASSERT_CALL(
    sleqp_sparse_vector_create_full(&fixed_initial,
                                    rosenbrock_num_variables - num_fixed));

  const double fixed_value
    = sleqp_sparse_vector_value_at(rosenbrock_initial, fixed_indices[0]);

  fixed_values[0] = fixed_value;

  ASSERT_CALL(sleqp_sparse_vector_push(fixed_initial, 0, fixed_value));

  ASSERT_CALL(sleqp_fixed_var_dyn_func_create(&fixed_var_func,
                                              dyn_rosenbrock_func,
                                              num_fixed,
                                              fixed_indices,
                                              fixed_values));

  ASSERT_CALL(sleqp_dyn_func_set_accuracy(dyn_rosenbrock_func, accuracy));

  ASSERT_CALL(sleqp_dyn_func_set_accuracy(fixed_var_func, accuracy));
}

void
teardown()
{
  ASSERT_CALL(sleqp_sparse_vector_free(&fixed_initial));

  ASSERT_CALL(sleqp_func_release(&fixed_var_func));

  dyn_rosenbrock_teardown();
}

START_TEST(test_func_val)
{
  bool reject;

  int func_grad_nnz, cons_val_nnz, cons_jac_nnz;

  ASSERT_CALL(sleqp_func_set_value(dyn_rosenbrock_func,
                                   rosenbrock_initial,
                                   SLEQP_VALUE_REASON_NONE,
                                   &reject,
                                   &func_grad_nnz,
                                   &cons_val_nnz,
                                   &cons_jac_nnz));

  double func_val;

  ASSERT_CALL(sleqp_func_val(dyn_rosenbrock_func, &func_val));

  assert(!reject);

  ASSERT_CALL(sleqp_func_set_value(fixed_var_func,
                                   fixed_initial,
                                   SLEQP_VALUE_REASON_NONE,
                                   &reject,
                                   &func_grad_nnz,
                                   &cons_val_nnz,
                                   &cons_jac_nnz));

  double fixed_func_val;

  ASSERT_CALL(sleqp_func_val(fixed_var_func, &fixed_func_val));

  ck_assert(sleqp_is_eq(func_val, fixed_func_val, 2. * accuracy));
}
END_TEST

START_TEST(test_func_grad)
{
}
END_TEST

Suite*
fixed_var_dynamic_test_suite()
{
  Suite* suite;
  TCase* tc_eval;

  suite = suite_create("Fixed variable dynamic function tests");

  tc_eval = tcase_create("Evaluation tests");

  tcase_add_checked_fixture(tc_eval, setup, teardown);

  tcase_add_test(tc_eval, test_func_val);

  tcase_add_test(tc_eval, test_func_grad);

  suite_add_tcase(suite, tc_eval);

  return suite;
}

TEST_MAIN(fixed_var_dynamic_test_suite)
