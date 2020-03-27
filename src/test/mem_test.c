#include <stdlib.h>
#include <limits.h>

#include <check.h>

#include "sleqp.h"

#include "sleqp.h"

#include "test_common.h"

START_TEST(test_alloc)
{
  int* ptr;

  ASSERT_CALL(sleqp_malloc(&ptr));

  ck_assert(ptr != NULL);

  sleqp_free(&ptr);

  ck_assert(ptr == NULL);
}
END_TEST

START_TEST(test_calloc)
{
  int* ptr = (int*) 1;

  ASSERT_CALL(sleqp_calloc(&ptr, 100));

  ck_assert(ptr != NULL);

  sleqp_free(&ptr);

  ck_assert(ptr == NULL);
}
END_TEST

START_TEST(test_calloc_zero)
{
  int* ptr = (int*) 1;

  ASSERT_CALL(sleqp_calloc(&ptr, 0));

  ck_assert(ptr == NULL);
}
END_TEST

/*
START_TEST(test_alloc_nomem)
{
  int* ptr = (int*) 1;

  SLEQP_RETCODE retcode = (sleqp_calloc(&ptr, INT_MAX / sizeof(int)));

  ck_assert_int_eq(retcode, SLEQP_NOMEM);

  ck_assert(ptr == NULL);
}
END_TEST
*/

START_TEST(test_realloc)
{
  int* ptr = (int*) 1;

  ASSERT_CALL(sleqp_calloc(&ptr, 100));

  ck_assert(ptr != NULL);

  ASSERT_CALL(sleqp_realloc(&ptr, 200));

  ck_assert(ptr != NULL);

  sleqp_free(&ptr);

  ck_assert(ptr == NULL);
}
END_TEST

/*
START_TEST(test_realloc_nomem)
{
  int* ptr = (int*) 1;

  ASSERT_CALL(sleqp_calloc(&ptr, 100));

  ck_assert(ptr != NULL);

  SLEQP_RETCODE retcode = (sleqp_realloc(&ptr, INT_MAX / sizeof(int)));

  ck_assert(ptr == NULL);

  ck_assert_int_eq(retcode, SLEQP_NOMEM);
}
END_TEST
*/

Suite* mem_test_suite()
{
  Suite *suite;
  TCase *tc_alloc;
  TCase *tc_realloc;

  suite = suite_create("Memory tests");

  tc_alloc = tcase_create("Allocation");

  tcase_add_test(tc_alloc, test_alloc);

  tcase_add_test(tc_alloc, test_calloc);

  tcase_add_test(tc_alloc, test_calloc_zero);

  // tcase_add_test(tc_alloc, test_alloc_nomem);

  tc_realloc = tcase_create("Reallocation");

  tcase_add_test(tc_realloc, test_realloc);

  // tcase_add_test(tc_alloc, test_realloc_nomem);

  suite_add_tcase(suite, tc_alloc);

  suite_add_tcase(suite, tc_realloc);

  return suite;
}

int main()
{
  int num_fails;
  Suite* suite;
  SRunner* srunner;

  suite =  mem_test_suite();
  srunner = srunner_create(suite);

  srunner_set_fork_status(srunner, CK_NOFORK);
  srunner_run_all(srunner, CK_NORMAL);

  num_fails = srunner_ntests_failed(srunner);

  srunner_free(srunner);

  return (num_fails > 0) ? EXIT_FAILURE : EXIT_SUCCESS;
}
