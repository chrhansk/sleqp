#include <check.h>
#include <stdlib.h>

#include "pub_settings.h"
#include "pub_types.h"
#include "settings.h"
#include "test_common.h"

const char* settings_file_name = "sleqp_example.opt";

START_TEST(test_read_settings)
{
  SleqpSettings* settings;

  ASSERT_CALL(sleqp_settings_create(&settings));

  ASSERT_CALL(sleqp_settings_read_file(settings, settings_file_name));

  ck_assert_double_eq(sleqp_settings_real_value(settings, SLEQP_SETTINGS_REAL_ZERO_EPS),
                      1e-2);

  ck_assert_int_eq(sleqp_settings_enum_value(settings, SLEQP_SETTINGS_ENUM_DUAL_ESTIMATION_TYPE),
                   SLEQP_DUAL_ESTIMATION_TYPE_LP);

  ck_assert_int_eq(sleqp_settings_int_value(settings, SLEQP_SETTINGS_INT_MAX_NEWTON_ITERATIONS),
                   10);

  ck_assert(sleqp_settings_bool_value(settings, SLEQP_SETTINGS_BOOL_GLOBAL_PENALTY_RESETS) == false);

  ASSERT_CALL(sleqp_settings_release(&settings));
}
END_TEST

Suite*
settings_test_suite()
{
  Suite* suite;
  TCase* tc_read;

  suite = suite_create("Settings tests");

  tc_read = tcase_create("Read settings");

  suite_add_tcase(suite, tc_read);

  tcase_add_test(tc_read, test_read_settings);

  return suite;
}

TEST_MAIN(settings_test_suite)
