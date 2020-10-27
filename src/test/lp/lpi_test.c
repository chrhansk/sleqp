#include <stdlib.h>
#include <check.h>

#include "sleqp_cmp.h"

#include "test_common.h"

#include "lp/sleqp_lpi.h"

START_TEST(test_simplex_solve)
{
  SleqpParams* params;

  SleqpLPi* lp_interface;

  SleqpSparseMatrix* cons_matrix;

  const double inf = sleqp_infinity();

  ASSERT_CALL(sleqp_params_create(&params));

  int num_variables = 2;
  int num_constraints = 1;

  ASSERT_CALL(sleqp_lpi_create_default_interface(&lp_interface,
                                                 num_variables,
                                                 num_constraints,
                                                 params));

  double objective[] = {-1, 0};
  double vars_lb[] = {0, 0};
  double vars_ub[] = {inf, inf};
  double cons_lb[] = {-inf};
  double cons_ub[] = {1};

  ASSERT_CALL(sleqp_sparse_matrix_create(&cons_matrix, 1, 2, 2));

  double* cons_data = sleqp_sparse_matrix_get_data(cons_matrix);
  int* cons_rows = sleqp_sparse_matrix_get_rows(cons_matrix);
  int* cons_cols = sleqp_sparse_matrix_get_cols(cons_matrix);

  cons_data[0] = 1;
  cons_data[1] = 1;

  cons_rows[0] = 0;
  cons_rows[1] = 0;

  cons_cols[0] = 0;
  cons_cols[1] = 1;
  cons_cols[2] = 2;

  ASSERT_CALL(sleqp_sparse_matrix_set_nnz(cons_matrix, 2));

  ASSERT_CALL(sleqp_lpi_set_bounds(lp_interface,
                                   cons_lb,
                                   cons_ub,
                                   vars_lb,
                                   vars_ub));

  ASSERT_CALL(sleqp_lpi_set_coefficients(lp_interface,
                                         cons_matrix));

  ASSERT_CALL(sleqp_lpi_set_objective(lp_interface,
                                      objective));

  ASSERT_CALL(sleqp_lpi_solve(lp_interface));

  double solution[] = {-1, -1};
  double objective_value = 0;

  ASSERT_CALL(sleqp_lpi_get_primal_sol(lp_interface,
                                       &objective_value,
                                       solution));

  double tolerance = 1e-8;

  ck_assert(sleqp_is_eq(objective_value, -1., tolerance));

  ck_assert(sleqp_is_eq(solution[0], 1., tolerance));
  ck_assert(sleqp_is_eq(solution[1], 0., tolerance));

  SLEQP_BASESTAT variable_stats[] = {0, 0};

  ASSERT_CALL(sleqp_lpi_get_varstats(lp_interface,
                                     variable_stats));

  ck_assert(variable_stats[0] == SLEQP_BASESTAT_BASIC);
  ck_assert(variable_stats[1] == SLEQP_BASESTAT_LOWER);

  double cons_dual = inf;
  double vars_dual[] = {inf, inf};

  ASSERT_CALL(sleqp_lpi_get_dual_sol(lp_interface,
                                     vars_dual,
                                     &cons_dual));

  ck_assert(sleqp_is_eq(cons_dual, -1., tolerance));

  ck_assert(sleqp_is_eq(vars_dual[0], 0., tolerance));
  ck_assert(sleqp_is_eq(vars_dual[1], 1., tolerance));

  ASSERT_CALL(sleqp_sparse_matrix_release(&cons_matrix));

  ASSERT_CALL(sleqp_lpi_free(&lp_interface));

  ASSERT_CALL(sleqp_params_free(&params));
}
END_TEST

Suite* lpi_suite()
{
  Suite *suite;
  TCase *tc_simplex;

  suite = suite_create("LPI");

  tc_simplex = tcase_create("Simplex solve");
  tcase_add_test(tc_simplex, test_simplex_solve);
  suite_add_tcase(suite, tc_simplex);

  return suite;
}

int main()
{
  int num_fails;
  Suite* suite;
  SRunner* srunner;

  suite = lpi_suite();
  srunner = srunner_create(suite);

  srunner_set_fork_status(srunner, CK_NOFORK);
  srunner_run_all(srunner, CK_NORMAL);

  num_fails = srunner_ntests_failed(srunner);

  srunner_free(srunner);

  return (num_fails > 0) ? EXIT_FAILURE : EXIT_SUCCESS;
}
