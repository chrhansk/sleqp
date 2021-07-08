#include <stdlib.h>
#include <check.h>

#include "test_common.h"
#include "quadcons_fixture.h"

#include "mem.h"
#include "util.h"
#include "sparse/sparse_matrix.h"

#include "preprocessor/fixed_var_func.h"

int num_fixed;
int* fixed_indices;
double* fixed_values;

static int num_variables;
static int num_constraints;

SleqpFunc* fixed_var_func;

SleqpSparseVec* value;
SleqpSparseVec* grad;

SleqpSparseVec* value;
SleqpSparseVec* fixed_value;

SleqpSparseVec* grad;
SleqpSparseVec* fixed_grad;

SleqpSparseMatrix* cons_jac;
SleqpSparseMatrix* fixed_cons_jac;

SleqpSparseVec* direction;
SleqpSparseVec* fixed_direction;

SleqpSparseVec* product;
SleqpSparseVec* fixed_product;

SleqpSparseVec* cons_duals;

void setup()
{
  quadconsfunc_setup();

  num_variables = quadconsfunc_num_variables;
  num_constraints = quadconsfunc_num_constraints;

  num_fixed = 1;

  ASSERT_CALL(sleqp_alloc_array(&fixed_indices, num_fixed));
  ASSERT_CALL(sleqp_alloc_array(&fixed_values, num_fixed));

  fixed_indices[0] = 0;
  fixed_values[0] = 1.;

  ASSERT_CALL(sleqp_fixed_var_func_create(&fixed_var_func,
                                          quadconsfunc,
                                          num_fixed,
                                          fixed_indices,
                                          fixed_values));

  ASSERT_CALL(sleqp_sparse_vector_create_full(&value,
                                              num_variables));

  ASSERT_CALL(sleqp_sparse_vector_create_full(&fixed_value,
                                              num_variables - num_fixed));

  ASSERT_CALL(sleqp_sparse_vector_create_full(&grad,
                                              num_variables));

  ASSERT_CALL(sleqp_sparse_vector_create_full(&fixed_grad,
                                              num_variables - num_fixed));

  ASSERT_CALL(sleqp_sparse_vector_push(fixed_value, 0, 2.));

  ASSERT_CALL(sleqp_sparse_vector_push(value, 0, 1.));
  ASSERT_CALL(sleqp_sparse_vector_push(value, 1, 2.));

  ASSERT_CALL(sleqp_sparse_matrix_create(&cons_jac,
                                         num_constraints,
                                         num_variables,
                                         0));

  ASSERT_CALL(sleqp_sparse_matrix_create(&fixed_cons_jac,
                                         num_constraints,
                                         num_variables - num_fixed,
                                         0));

  ASSERT_CALL(sleqp_sparse_vector_create_full(&direction,
                                              num_variables));

  // This value should be stripped during evaluation
  ASSERT_CALL(sleqp_sparse_vector_push(direction, 0, 10.));
  ASSERT_CALL(sleqp_sparse_vector_push(direction, 1, 1.));

  ASSERT_CALL(sleqp_sparse_vector_create_full(&fixed_direction,
                                              num_variables - num_fixed));

  ASSERT_CALL(sleqp_sparse_vector_push(fixed_direction, 0, 1.));

  ASSERT_CALL(sleqp_sparse_vector_create_full(&product,
                                              num_variables));

  ASSERT_CALL(sleqp_sparse_vector_create_full(&fixed_product,
                                              num_variables - num_fixed));

  ASSERT_CALL(sleqp_sparse_vector_create_full(&cons_duals,
                                              num_constraints));

  for(int i = 0; i < num_constraints; ++i)
  {
    ASSERT_CALL(sleqp_sparse_vector_push(cons_duals, i, 1.));
  }
}

void teardown()
{
  ASSERT_CALL(sleqp_sparse_vector_free(&cons_duals));

  ASSERT_CALL(sleqp_sparse_vector_free(&fixed_product));
  ASSERT_CALL(sleqp_sparse_vector_free(&product));

  ASSERT_CALL(sleqp_sparse_vector_free(&fixed_direction));
  ASSERT_CALL(sleqp_sparse_vector_free(&direction));

  ASSERT_CALL(sleqp_sparse_matrix_release(&fixed_cons_jac));
  ASSERT_CALL(sleqp_sparse_matrix_release(&cons_jac));

  ASSERT_CALL(sleqp_sparse_vector_free(&fixed_grad));
  ASSERT_CALL(sleqp_sparse_vector_free(&grad));

  ASSERT_CALL(sleqp_sparse_vector_free(&fixed_value));
  ASSERT_CALL(sleqp_sparse_vector_free(&value));

  ASSERT_CALL(sleqp_func_release(&fixed_var_func));

  sleqp_free(&fixed_values);
  sleqp_free(&fixed_indices);

  quadconsfunc_teardown();
}

START_TEST(test_func_eval)
{
  int func_grad_nnz, cons_val_nnz, cons_jac_nnz;

  ASSERT_CALL(sleqp_func_set_value(fixed_var_func,
                                   fixed_value,
                                   SLEQP_VALUE_REASON_NONE,
                                   &func_grad_nnz,
                                   &cons_val_nnz,
                                   &cons_jac_nnz));

  double fixed_func_val;

  ASSERT_CALL(sleqp_func_val(fixed_var_func, &fixed_func_val));

  ASSERT_CALL(sleqp_func_set_value(quadconsfunc,
                                   value,
                                   SLEQP_VALUE_REASON_NONE,
                                   &func_grad_nnz,
                                   &cons_val_nnz,
                                   &cons_jac_nnz));

  double func_val;

  ASSERT_CALL(sleqp_func_val(quadconsfunc, &func_val));

  ck_assert(fixed_func_val == func_val);

}
END_TEST

START_TEST(test_func_grad)
{
  int func_grad_nnz, cons_val_nnz, cons_jac_nnz;

  ASSERT_CALL(sleqp_func_set_value(fixed_var_func,
                                   fixed_value,
                                   SLEQP_VALUE_REASON_NONE,
                                   &func_grad_nnz,
                                   &cons_val_nnz,
                                   &cons_jac_nnz));

  ASSERT_CALL(sleqp_func_grad(fixed_var_func, fixed_grad));

  ASSERT_CALL(sleqp_func_set_value(quadconsfunc,
                                   value,
                                   SLEQP_VALUE_REASON_NONE,
                                   &func_grad_nnz,
                                   &cons_val_nnz,
                                   &cons_jac_nnz));

  ASSERT_CALL(sleqp_func_grad(quadconsfunc, grad));

  ck_assert(sleqp_sparse_vector_value_at(fixed_grad, 0) ==
            sleqp_sparse_vector_value_at(grad, 1));

}
END_TEST

START_TEST(test_cons_jac)
{
  int func_grad_nnz, cons_val_nnz, cons_jac_nnz;

  ASSERT_CALL(sleqp_func_set_value(fixed_var_func,
                                   fixed_value,
                                   SLEQP_VALUE_REASON_NONE,
                                   &func_grad_nnz,
                                   &cons_val_nnz,
                                   &cons_jac_nnz));

  ASSERT_CALL(sleqp_sparse_matrix_reserve(cons_jac,
                                          cons_jac_nnz));

  ASSERT_CALL(sleqp_func_cons_jac(fixed_var_func,
                                  NULL,
                                  fixed_cons_jac));

  ASSERT_CALL(sleqp_func_set_value(quadconsfunc,
                                   value,
                                   SLEQP_VALUE_REASON_NONE,
                                   &func_grad_nnz,
                                   &cons_val_nnz,
                                   &cons_jac_nnz));

  ASSERT_CALL(sleqp_sparse_matrix_reserve(cons_jac,
                                          cons_jac_nnz));

  ASSERT_CALL(sleqp_func_cons_jac(quadconsfunc,
                                  NULL,
                                  cons_jac));

  for(int i = 0; i < num_constraints; ++i)
  {
    ck_assert(sleqp_sparse_matrix_value_at(cons_jac, i, 1) ==
              sleqp_sparse_matrix_value_at(fixed_cons_jac, i, 0));
  }
}
END_TEST

START_TEST(test_hess_prod)
{
  int func_grad_nnz, cons_val_nnz, cons_jac_nnz;

  ASSERT_CALL(sleqp_func_set_value(fixed_var_func,
                                   fixed_value,
                                   SLEQP_VALUE_REASON_NONE,
                                   &func_grad_nnz,
                                   &cons_val_nnz,
                                   &cons_jac_nnz));

  const double one = 1.;

  ASSERT_CALL(sleqp_func_hess_prod(fixed_var_func,
                                   &one,
                                   fixed_direction,
                                   cons_duals,
                                   fixed_product));

  ASSERT_CALL(sleqp_func_set_value(quadconsfunc,
                                   value,
                                   SLEQP_VALUE_REASON_NONE,
                                   &func_grad_nnz,
                                   &cons_val_nnz,
                                   &cons_jac_nnz));

  ASSERT_CALL(sleqp_func_hess_prod(quadconsfunc,
                                   &one,
                                   direction,
                                   cons_duals,
                                   product));

  ck_assert(sleqp_sparse_vector_value_at(fixed_product, 0) ==
            sleqp_sparse_vector_value_at(product, 1));

}
END_TEST

Suite* fixed_var_test_suite()
{
  Suite *suite;
  TCase *tc_eval;

  suite = suite_create("Fixed variable function tests");

  tc_eval = tcase_create("Evaluation tests");

  tcase_add_checked_fixture(tc_eval,
                            setup,
                            teardown);

  tcase_add_test(tc_eval, test_func_eval);

  tcase_add_test(tc_eval, test_func_grad);

  tcase_add_test(tc_eval, test_cons_jac);

  tcase_add_test(tc_eval, test_hess_prod);

  suite_add_tcase(suite, tc_eval);

  return suite;
}

TEST_MAIN(fixed_var_test_suite)
