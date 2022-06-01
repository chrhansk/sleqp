#include <check.h>
#include <stdlib.h>

#include "cmp.h"

#include "test_common.h"

#include "lp/lpi.h"

SleqpParams* params;
SleqpOptions* options;

const int num_variables   = 2;
const int num_constraints = 1;

SleqpLPi* lp_interface;

SleqpSparseMatrix* cons_matrix;

void
set_lp_data()
{
  const double inf = sleqp_infinity();

  double objective[] = {-1, 0};
  double vars_lb[]   = {0, 0};
  double vars_ub[]   = {inf, inf};
  double cons_lb[]   = {-inf};
  double cons_ub[]   = {1};

  ASSERT_CALL(
    sleqp_lpi_set_bounds(lp_interface, cons_lb, cons_ub, vars_lb, vars_ub));

  ASSERT_CALL(sleqp_lpi_set_coeffs(lp_interface, cons_matrix));

  ASSERT_CALL(sleqp_lpi_set_objective(lp_interface, objective));
}

void
lpi_setup()
{
  ASSERT_CALL(sleqp_params_create(&params));
  ASSERT_CALL(sleqp_options_create(&options));

  ASSERT_CALL(sleqp_lpi_create_default(&lp_interface,
                                       num_variables,
                                       num_constraints,
                                       params,
                                       options));

  ASSERT_CALL(sleqp_sparse_matrix_create(&cons_matrix, 1, 2, 2));

  double* cons_data = sleqp_sparse_matrix_data(cons_matrix);
  int* cons_rows    = sleqp_sparse_matrix_rows(cons_matrix);
  int* cons_cols    = sleqp_sparse_matrix_cols(cons_matrix);

  cons_data[0] = 1;
  cons_data[1] = 1;

  cons_rows[0] = 0;
  cons_rows[1] = 0;

  cons_cols[0] = 0;
  cons_cols[1] = 1;
  cons_cols[2] = 2;

  ASSERT_CALL(sleqp_sparse_matrix_set_nnz(cons_matrix, 2));

  set_lp_data();
}

void
lpi_teardown()
{
  ASSERT_CALL(sleqp_sparse_matrix_release(&cons_matrix));

  ASSERT_CALL(sleqp_lpi_release(&lp_interface));

  ASSERT_CALL(sleqp_options_release(&options));
  ASSERT_CALL(sleqp_params_release(&params));
}

START_TEST(test_solve)
{
  ASSERT_CALL(sleqp_lpi_solve(lp_interface));

  double solution[]      = {-1, -1};
  double objective_value = 0;

  ASSERT_CALL(sleqp_lpi_primal_sol(lp_interface, &objective_value, solution));

  double tolerance = 1e-8;

  ck_assert(sleqp_is_eq(objective_value, -1., tolerance));

  ck_assert(sleqp_is_eq(solution[0], 1., tolerance));
  ck_assert(sleqp_is_eq(solution[1], 0., tolerance));

  SLEQP_BASESTAT var_stats[] = {0, 0};

  ASSERT_CALL(sleqp_lpi_vars_stats(lp_interface, var_stats));

  ck_assert(var_stats[0] == SLEQP_BASESTAT_BASIC);
  ck_assert(var_stats[1] == SLEQP_BASESTAT_LOWER);

  double cons_dual;
  double vars_dual[2];

  ASSERT_CALL(sleqp_lpi_dual_sol(lp_interface, vars_dual, &cons_dual));

  ck_assert(sleqp_is_eq(cons_dual, -1., tolerance));

  ck_assert(sleqp_is_eq(vars_dual[0], 0., tolerance));
  ck_assert(sleqp_is_eq(vars_dual[1], 1., tolerance));
}
END_TEST

START_TEST(test_basis_roundtrip)
{
  ASSERT_CALL(sleqp_lpi_solve(lp_interface));

  SLEQP_BASESTAT var_stats[num_variables];
  SLEQP_BASESTAT cons_stats[num_constraints];

  ASSERT_CALL(sleqp_lpi_vars_stats(lp_interface, var_stats));

  ASSERT_CALL(sleqp_lpi_cons_stats(lp_interface, cons_stats));

  set_lp_data();

  ASSERT_CALL(sleqp_lpi_set_basis(lp_interface, 0, var_stats, cons_stats));

  ASSERT_CALL(sleqp_lpi_restore_basis(lp_interface, 0));

  ASSERT_CALL(sleqp_lpi_solve(lp_interface));
}
END_TEST

Suite*
lpi_suite()
{
  Suite* suite;
  TCase* tc_solve;

  suite = suite_create("LP interface tests");

  tc_solve = tcase_create("LP interface solution");

  tcase_add_checked_fixture(tc_solve, lpi_setup, lpi_teardown);

  // tcase_add_test(tc_solve, test_solve);
  tcase_add_test(tc_solve, test_basis_roundtrip);

  suite_add_tcase(suite, tc_solve);

  return suite;
}

TEST_MAIN(lpi_suite)
