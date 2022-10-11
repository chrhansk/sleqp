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

SleqpVec* fixed_initial;

void
setup()
{
  dyn_rosenbrock_setup();

  ASSERT_CALL(
    sleqp_vec_create_full(&fixed_initial, rosenbrock_num_vars - num_fixed));

  const double fixed_value
    = sleqp_vec_value_at(rosenbrock_initial, fixed_indices[0]);

  fixed_values[0] = fixed_value;

  ASSERT_CALL(sleqp_vec_push(fixed_initial, 0, fixed_value));

  ASSERT_CALL(sleqp_fixed_var_dyn_func_create(&fixed_var_func,
                                              dyn_rosenbrock_func,
                                              num_fixed,
                                              fixed_indices,
                                              fixed_values));

  ASSERT_CALL(sleqp_dyn_func_set_error_bound(dyn_rosenbrock_func, accuracy));
  ASSERT_CALL(sleqp_dyn_func_set_obj_weight(dyn_rosenbrock_func, 1.));

  ASSERT_CALL(sleqp_dyn_func_set_error_bound(fixed_var_func, accuracy));
  ASSERT_CALL(sleqp_dyn_func_set_obj_weight(fixed_var_func, 1.));
}

void
teardown()
{
  ASSERT_CALL(sleqp_vec_free(&fixed_initial));

  ASSERT_CALL(sleqp_func_release(&fixed_var_func));

  dyn_rosenbrock_teardown();
}

START_TEST(test_obj_val)
{
  bool reject;

  ASSERT_CALL(sleqp_func_set_value(dyn_rosenbrock_func,
                                   rosenbrock_initial,
                                   SLEQP_VALUE_REASON_NONE,
                                   &reject));

  double obj_val;

  ASSERT_CALL(sleqp_func_obj_val(dyn_rosenbrock_func, &obj_val));

  ck_assert(!reject);

  ASSERT_CALL(sleqp_func_set_value(fixed_var_func,
                                   fixed_initial,
                                   SLEQP_VALUE_REASON_NONE,
                                   &reject));

  ck_assert(!reject);

  double fixed_obj_val;

  ASSERT_CALL(sleqp_func_obj_val(fixed_var_func, &fixed_obj_val));

  ck_assert(sleqp_is_eq(obj_val, fixed_obj_val, 2. * accuracy));
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

  tcase_add_test(tc_eval, test_obj_val);

  suite_add_tcase(suite, tc_eval);

  return suite;
}

TEST_MAIN(fixed_var_dynamic_test_suite)
