#include <stdlib.h>
#include <check.h>

#include "test_common.h"

#include "sparse/sleqp_sparse.h"


START_TEST(test_sparse_construction)
{
  SleqpSparseMatrix* identity;

  size_t size = 5;

  ASSERT_CALL(sleqp_sparse_matrix_create(&identity,
                                         5,
                                         5,
                                         5));

  for(size_t current = 0; current < size; ++current)
  {
    ASSERT_CALL(sleqp_sparse_matrix_add_column(identity, current));

    ASSERT_CALL(sleqp_sparse_matrix_push(identity,
                                         current,
                                         current,
                                         1.));
  }

  for(size_t row = 0; row < identity->num_rows; ++row)
  {
    for(size_t col = 0; col < identity->num_cols; ++col)
    {
      double* value = sleqp_sparse_matrix_at(identity, row, col);

      if(row != col)
      {
        ck_assert(!value);
      }
      else
      {
        ck_assert(value);
        ck_assert(*value == 1.);
      }
    }
  }


  ASSERT_CALL(sleqp_sparse_matrix_free(&identity));

}
END_TEST


Suite* sparse_test_suite()
{
  Suite *suite;
  TCase *tc_sparse_construction;

  suite = suite_create("Sparse tests");

  tc_sparse_construction = tcase_create("Sparse matrix construction");

  tcase_add_test(tc_sparse_construction, test_sparse_construction);
  suite_add_tcase(suite, tc_sparse_construction);

  return suite;
}


int main()
{
  int num_fails;
  Suite* suite;
  SRunner* srunner;

  suite = sparse_test_suite();
  srunner = srunner_create(suite);

  srunner_run_all(srunner, CK_NORMAL);

  num_fails = srunner_ntests_failed(srunner);

  srunner_free(srunner);

  return (num_fails > 0) ? EXIT_FAILURE : EXIT_SUCCESS;
}
