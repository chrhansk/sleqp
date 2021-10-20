#include <check.h>
#include <string.h>

#include "test_common.h"

bool handler_called = false;

const char* handler_message;

const char* error_message = "error";
const char* info_message = "info";

void log_handler(SLEQP_LOG_LEVEL level,
                 time_t time,
                 const char* message)
{
  handler_message = message;
  handler_called = true;
}

void log_setup()
{
  handler_called = false;

  sleqp_log_set_handler(log_handler);
}

void log_teardown()
{

}

START_TEST(test_log_msg)
{
  sleqp_log_set_level(SLEQP_LOG_ERROR);

  sleqp_log_error(error_message);
  sleqp_log_info(info_message);

  ck_assert(handler_called);

  ck_assert(!strcmp(handler_message, error_message));
}
END_TEST

START_TEST(test_log_trace)
{
  sleqp_log_set_level(SLEQP_LOG_ERROR);

  sleqp_log_trace_level(SLEQP_LOG_ERROR, "<file>", 1, "<message>");

  ck_assert(handler_called);

  ck_assert(strstr(handler_message, "<file>"));
  ck_assert(strstr(handler_message, "1"));
  ck_assert(strstr(handler_message, "<message>"));
}
END_TEST

START_TEST(test_silent)
{
  sleqp_log_set_level(SLEQP_LOG_SILENT);

  sleqp_log_error(error_message);
  sleqp_log_info(info_message);

  ck_assert(!handler_called);
}
END_TEST

Suite* log_test_suite()
{
  Suite *suite;
  TCase *tc_log;

  suite = suite_create("Logging tests");

  tc_log = tcase_create("Logging test");

  tcase_add_checked_fixture(tc_log,
                            log_setup,
                            log_teardown);

  tcase_add_test(tc_log, test_silent);

  tcase_add_test(tc_log, test_log_msg);

  tcase_add_test(tc_log, test_log_trace);

  suite_add_tcase(suite, tc_log);

  return suite;
}

TEST_MAIN(log_test_suite)
