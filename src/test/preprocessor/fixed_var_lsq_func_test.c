#include <check.h>

#include "test_common.h"

#include "rosenbrock_fixture.h"
#include "rosenbrock_lsq_fixture.h"

#include "lsq.h"
#include "params.h"

#include "preprocessor/fixed_var_func.h"

SleqpSparseVec* residuals;
SleqpSparseVec* fixed_residuals;

SleqpSparseVec* forward;
SleqpSparseVec* fixed_forward;

SleqpSparseVec* adjoint;
SleqpSparseVec* fixed_adjoint;

SleqpParams* params;

int num_fixed = 1;
int fixed_indices[] = {0};
double fixed_values[] = {1.};

SleqpFunc* fixed_var_func;

SleqpSparseVec* fixed_initial;

void setup()
{
  rosenbrock_setup();

  rosenbrock_lsq_setup();

  ASSERT_CALL(sleqp_params_create(&params));

  ASSERT_CALL(sleqp_sparse_vector_create_empty(&residuals,
                                               rosenbrock_num_residuals));

  ASSERT_CALL(sleqp_sparse_vector_create_empty(&fixed_residuals,
                                               rosenbrock_num_residuals));

  ASSERT_CALL(sleqp_sparse_vector_create_full(&forward,
                                              rosenbrock_num_variables));

  ASSERT_CALL(sleqp_sparse_vector_create_full(&fixed_forward,
                                              rosenbrock_num_variables - num_fixed));

  ASSERT_CALL(sleqp_sparse_vector_create_full(&adjoint,
                                              rosenbrock_num_residuals));

  ASSERT_CALL(sleqp_sparse_vector_create_full(&fixed_adjoint,
                                              rosenbrock_num_residuals));

  fixed_values[0] = sleqp_sparse_vector_value_at(rosenbrock_initial,
                                                 fixed_indices[0]);

  ASSERT_CALL(sleqp_fixed_var_lsq_func_create(&fixed_var_func,
                                              rosenbrock_lsq_func,
                                              params,
                                              num_fixed,
                                              fixed_indices,
                                              fixed_values));

  ASSERT_CALL(sleqp_sparse_vector_create_full(&fixed_initial,
                                              rosenbrock_num_variables - num_fixed));

  const double fixed_value = sleqp_sparse_vector_value_at(rosenbrock_initial,
                                                          1);

  fixed_values[0] = fixed_value;

  ASSERT_CALL(sleqp_sparse_vector_push(fixed_initial,
                                       0,
                                       fixed_value));

  bool reject;

  int func_grad_nnz, cons_val_nnz, cons_jac_nnz;

  ASSERT_CALL(sleqp_func_set_value(rosenbrock_lsq_func,
                                   rosenbrock_initial,
                                   SLEQP_VALUE_REASON_NONE,
                                   &reject,
                                   &func_grad_nnz,
                                   &cons_val_nnz,
                                   &cons_jac_nnz));

  assert(!reject);

  ASSERT_CALL(sleqp_func_set_value(fixed_var_func,
                                   fixed_initial,
                                   SLEQP_VALUE_REASON_NONE,
                                   &reject,
                                   &func_grad_nnz,
                                   &cons_val_nnz,
                                   &cons_jac_nnz));

  assert(!reject);
}

void teardown()
{
  ASSERT_CALL(sleqp_sparse_vector_free(&fixed_initial));

  ASSERT_CALL(sleqp_func_release(&fixed_var_func));

  ASSERT_CALL(sleqp_sparse_vector_free(&fixed_adjoint));
  ASSERT_CALL(sleqp_sparse_vector_free(&adjoint));

  ASSERT_CALL(sleqp_sparse_vector_free(&fixed_forward));
  ASSERT_CALL(sleqp_sparse_vector_free(&forward));

  ASSERT_CALL(sleqp_sparse_vector_free(&fixed_residuals));
  ASSERT_CALL(sleqp_sparse_vector_free(&residuals));

  ASSERT_CALL(sleqp_params_release(&params));

  rosenbrock_lsq_teardown();

  rosenbrock_teardown();
}

START_TEST(test_residuals)
{
  const double zero_eps = sleqp_params_get(params, SLEQP_PARAM_ZERO_EPS);

  ASSERT_CALL(sleqp_lsq_func_residuals(rosenbrock_lsq_func,
                                       residuals));

  ASSERT_CALL(sleqp_lsq_func_residuals(fixed_var_func,
                                       fixed_residuals));

  ck_assert(sleqp_sparse_vector_eq(residuals,
                                   fixed_residuals,
                                   zero_eps));
}
END_TEST

START_TEST(test_jac_forward)
{
  const double zero_eps = sleqp_params_get(params, SLEQP_PARAM_ZERO_EPS);

  ASSERT_CALL(sleqp_sparse_vector_push(forward, 1, 1.));

  ASSERT_CALL(sleqp_lsq_func_jac_forward(rosenbrock_lsq_func,
                                         forward,
                                         adjoint));

  ASSERT_CALL(sleqp_sparse_vector_push(fixed_forward, 0, 1.));

  ASSERT_CALL(sleqp_lsq_func_jac_forward(fixed_var_func,
                                         fixed_forward,
                                         fixed_adjoint));

  ck_assert(sleqp_sparse_vector_eq(adjoint,
                                   fixed_adjoint,
                                   zero_eps));
}
END_TEST

START_TEST(test_jac_adjoint)
{
  ASSERT_CALL(sleqp_sparse_vector_push(adjoint, 1, 1.));

  ASSERT_CALL(sleqp_lsq_func_jac_adjoint(rosenbrock_lsq_func,
                                         adjoint,
                                         forward));

  ASSERT_CALL(sleqp_lsq_func_jac_adjoint(fixed_var_func,
                                         adjoint,
                                         fixed_forward));

  ck_assert(sleqp_sparse_vector_value_at(forward, 1) ==
            sleqp_sparse_vector_value_at(fixed_forward, 0));
}
END_TEST

Suite* fixed_var_lsq_test_suite()
{
  Suite *suite;
  TCase *tc_eval;

  suite = suite_create("Fixed variable LSQ function tests");

  tc_eval = tcase_create("Evaluation tests");

  tcase_add_checked_fixture(tc_eval,
                            setup,
                            teardown);

  tcase_add_test(tc_eval, test_residuals);

  tcase_add_test(tc_eval, test_jac_forward);

  tcase_add_test(tc_eval, test_jac_adjoint);

  suite_add_tcase(suite, tc_eval);

  return suite;
}

TEST_MAIN(fixed_var_lsq_test_suite)
