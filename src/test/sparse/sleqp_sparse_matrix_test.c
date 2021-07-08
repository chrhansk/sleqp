#include <stdlib.h>
#include <check.h>

#include "test_common.h"

#include "cmp.h"
#include "mem.h"

#include "sparse/sparse_matrix.h"
#include "sparse/sparse_vec.h"

START_TEST(test_sparse_matrix_vector_product)
{
  SleqpSparseMatrix* matrix;

  int num_rows = 2;
  int num_cols = 3;
  int nnz_max = 4;

  ASSERT_CALL(sleqp_sparse_matrix_create(&matrix,
                                         num_rows,
                                         num_cols,
                                         nnz_max));

  // row, col, val
  ASSERT_CALL(sleqp_sparse_matrix_push_column(matrix, 0));
  ASSERT_CALL(sleqp_sparse_matrix_push(matrix, 0, 0, 1.));

  ASSERT_CALL(sleqp_sparse_matrix_push_column(matrix, 1));
  ASSERT_CALL(sleqp_sparse_matrix_push(matrix, 1, 1, 2.));

  ASSERT_CALL(sleqp_sparse_matrix_push_column(matrix, 2));
  ASSERT_CALL(sleqp_sparse_matrix_push(matrix, 0, 2, 2.));
  ASSERT_CALL(sleqp_sparse_matrix_push(matrix, 1, 2, 3.));

  SleqpSparseVec* vector;
  double* result;

  ASSERT_CALL(sleqp_sparse_vector_create(&vector, 3, 3));
  ASSERT_CALL(sleqp_sparse_vector_push(vector, 0, 2.));
  ASSERT_CALL(sleqp_sparse_vector_push(vector, 1, 4.));
  ASSERT_CALL(sleqp_sparse_vector_push(vector, 2, 3.));

  ASSERT_CALL(sleqp_alloc_array(&result, 2));

  ASSERT_CALL(sleqp_sparse_matrix_vector_product(matrix,
                                                 vector,
                                                 result));

  double tolerance = 1e-8;

  ck_assert(sleqp_is_eq(result[0], 8., tolerance));
  ck_assert(sleqp_is_eq(result[1], 17., tolerance));

  sleqp_free(&result);

  ASSERT_CALL(sleqp_sparse_vector_free(&vector));

  ASSERT_CALL(sleqp_sparse_matrix_release(&matrix));

}
END_TEST

START_TEST(test_sparse_reserve)
{
  SleqpSparseMatrix* matrix;

  int size = 5;

  ASSERT_CALL(sleqp_sparse_matrix_create(&matrix,
                                         size,
                                         size,
                                         0));

  ASSERT_CALL(sleqp_sparse_matrix_reserve(matrix, 10));

  ck_assert(sleqp_sparse_matrix_get_nnz_max(matrix) >= 10);
  ck_assert_int_eq(sleqp_sparse_matrix_get_nnz(matrix), 0);

  ASSERT_CALL(sleqp_sparse_matrix_release(&matrix));
}
END_TEST

START_TEST(test_sparse_increase_size)
{
  SleqpSparseMatrix* matrix;

  int num_nnz = 5;
  int initial_size = 3, size = 10;

  ASSERT_CALL(sleqp_sparse_matrix_create(&matrix,
                                         initial_size,
                                         initial_size,
                                         0));

  ASSERT_CALL(sleqp_sparse_matrix_reserve(matrix, num_nnz));

  ASSERT_CALL(sleqp_sparse_matrix_resize(matrix, size, size));

  ck_assert_int_eq(sleqp_sparse_matrix_get_num_rows(matrix), size);
  ck_assert_int_eq(sleqp_sparse_matrix_get_num_cols(matrix), size);

  const int* matrix_cols = sleqp_sparse_matrix_get_cols(matrix);

  for(int col = 0; col < size; ++col)
  {
    ck_assert_int_eq(matrix_cols[col], 0);
  }

  ck_assert_int_eq(sleqp_sparse_matrix_get_nnz(matrix), 0);
  ck_assert_int_eq(sleqp_sparse_matrix_get_nnz_max(matrix), num_nnz);

  ASSERT_CALL(sleqp_sparse_matrix_release(&matrix));

}
END_TEST

START_TEST(test_sparse_pop_column)
{
  SleqpSparseMatrix* matrix;

  int size = 5;

  ASSERT_CALL(sleqp_sparse_matrix_create(&matrix,
                                         size,
                                         size,
                                         size));

  for(int current = 0; current < size; ++current)
  {
    ASSERT_CALL(sleqp_sparse_matrix_push_column(matrix, current));

    ASSERT_CALL(sleqp_sparse_matrix_push(matrix,
                                         current,
                                         current,
                                         1.));
  }

  ck_assert_int_eq(sleqp_sparse_matrix_get_nnz(matrix), size);

  int removed = 0;

  for(int column = size - 1; column > 0; --column)
  {
    ASSERT_CALL(sleqp_sparse_matrix_pop_column(matrix, column));

    ++removed;

    ck_assert_int_eq(sleqp_sparse_matrix_get_nnz(matrix), size - removed);

  }


  ASSERT_CALL(sleqp_sparse_matrix_release(&matrix));

}
END_TEST

START_TEST(test_sparse_construction)
{
  SleqpSparseMatrix* identity;

  int size = 5;

  ASSERT_CALL(sleqp_sparse_matrix_create(&identity,
                                         size,
                                         size,
                                         size));

  ck_assert_int_eq(sleqp_sparse_matrix_get_num_cols(identity), size);
  ck_assert_int_eq(sleqp_sparse_matrix_get_num_rows(identity), size);
  ck_assert_int_eq(sleqp_sparse_matrix_get_nnz_max(identity), size);
  ck_assert_int_eq(sleqp_sparse_matrix_get_nnz(identity), 0);

  for(int current = 0; current < size; ++current)
  {
    ASSERT_CALL(sleqp_sparse_matrix_push_column(identity, current));

    ASSERT_CALL(sleqp_sparse_matrix_push(identity,
                                         current,
                                         current,
                                         1.));

    ck_assert_int_eq(sleqp_sparse_matrix_get_nnz(identity), current + 1);
  }

  for(int row = 0; row < size; ++row)
  {
    for(int col = 0; col < size; ++col)
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

  ASSERT_CALL(sleqp_sparse_matrix_release(&identity));
}
END_TEST

START_TEST(test_sparse_decrease_size)
{
  SleqpSparseMatrix* identity;

  int size = 5, reduced_size = 2;

  ASSERT_CALL(sleqp_sparse_matrix_create(&identity,
                                         size,
                                         size,
                                         size));

  for(int current = 0; current < size; ++current)
  {
    ASSERT_CALL(sleqp_sparse_matrix_push_column(identity, current));

    ASSERT_CALL(sleqp_sparse_matrix_push(identity,
                                         current,
                                         current,
                                         1.));
  }

  ck_assert_int_eq(sleqp_sparse_matrix_get_nnz(identity), size);

  ASSERT_CALL(sleqp_sparse_matrix_resize(identity,
                                         reduced_size,
                                         reduced_size));

  ck_assert_int_eq(sleqp_sparse_matrix_get_nnz(identity), reduced_size);

  ASSERT_CALL(sleqp_sparse_matrix_release(&identity));
}
END_TEST

Suite* sparse_test_suite()
{
  Suite *suite;
  TCase *tc_sparse_construction;
  TCase *tc_sparse_modification;
  TCase *tc_sparse_operations;

  suite = suite_create("Sparse tests");

  tc_sparse_construction = tcase_create("Sparse matrix construction");

  tcase_add_test(tc_sparse_construction, test_sparse_construction);

  tc_sparse_modification = tcase_create("Sparse matrix modification");

  tc_sparse_operations = tcase_create("Sparse matrix operations");

  tcase_add_test(tc_sparse_modification, test_sparse_reserve);
  tcase_add_test(tc_sparse_modification, test_sparse_increase_size);
  tcase_add_test(tc_sparse_modification, test_sparse_decrease_size);
  tcase_add_test(tc_sparse_modification, test_sparse_pop_column);

  tcase_add_test(tc_sparse_operations, test_sparse_matrix_vector_product);

  suite_add_tcase(suite, tc_sparse_construction);

  suite_add_tcase(suite, tc_sparse_modification);

  suite_add_tcase(suite, tc_sparse_operations);

  return suite;
}

TEST_MAIN(sparse_test_suite)
