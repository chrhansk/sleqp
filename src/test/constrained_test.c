#include <stdlib.h>
#include <check.h>

#include "sleqp.h"
#include "sleqp_aug_jacobian.h"
#include "sleqp_cauchy.h"
#include "sleqp_cmp.h"
#include "sleqp_dual_estimation.h"
#include "sleqp_mem.h"

#include "lp/sleqp_lpi_soplex.h"

#include "test_common.h"

#include "quadcons_fixture.h"

SleqpProblem* problem;
SleqpIterate* iterate;
SleqpLPi* lp_interface;
SleqpCauchyData* cauchy_data;

SleqpSparseVec* cauchy_direction;

void constrained_setup()
{
  quadconsfunc_setup();

  ASSERT_CALL(sleqp_problem_create(&problem,
                                   quadconsfunc,
                                   quadconsfunc_var_lb,
                                   quadconsfunc_var_ub,
                                   quadconsfunc_cons_lb,
                                   quadconsfunc_cons_ub));

  ASSERT_CALL(sleqp_iterate_create(&iterate,
                                   problem,
                                   quadconsfunc_x));

  int num_lp_variables = problem->num_variables + 2*problem->num_constraints;
  int num_lp_constraints = problem->num_constraints;

  ASSERT_CALL(sleqp_lpi_soplex_create_interface(&lp_interface,
                                                num_lp_variables,
                                                num_lp_constraints));

  ASSERT_CALL(sleqp_set_and_evaluate(problem, iterate));

  ASSERT_CALL(sleqp_cauchy_data_create(&cauchy_data,
                                       problem,
                                       lp_interface));

  ASSERT_CALL(sleqp_sparse_vector_create(&cauchy_direction, 2, 2));

}

START_TEST(test_active_set)
{
  SleqpActiveSet* active_set = iterate->active_set;

  double trust_radius = 0.1, penalty_parameter = 1.;

  ASSERT_CALL(sleqp_cauchy_set_iterate(cauchy_data,
                                       iterate,
                                       trust_radius));

  ASSERT_CALL(sleqp_cauchy_solve(cauchy_data,
                                 iterate->func_grad,
                                 penalty_parameter));

  ASSERT_CALL(sleqp_cauchy_get_active_set(cauchy_data,
                                          iterate,
                                          trust_radius));

  SLEQP_ACTIVE_STATE* cons_states = sleqp_active_set_cons_states(active_set);

  ck_assert_int_eq(cons_states[0], SLEQP_INACTIVE);
  ck_assert_int_eq(cons_states[1], SLEQP_ACTIVE_UPPER);
}
END_TEST

START_TEST(test_dual_variable)
{
  SleqpSparseVec* cons_dual = iterate->cons_dual;

  SleqpAugJacobian* jacobian;
  SleqpDualEstimationData* estimation_data;

  double trust_radius = 0.1, penalty_parameter = 1.;

  ASSERT_CALL(sleqp_cauchy_set_iterate(cauchy_data,
                                       iterate,
                                       trust_radius));

  ASSERT_CALL(sleqp_cauchy_solve(cauchy_data,
                                 iterate->func_grad,
                                 penalty_parameter));

  ASSERT_CALL(sleqp_cauchy_get_active_set(cauchy_data,
                                          iterate,
                                          trust_radius));

  ASSERT_CALL(sleqp_cauchy_get_direction(cauchy_data, iterate, cauchy_direction));

  ASSERT_CALL(sleqp_aug_jacobian_create(&jacobian,
                                        problem));

  ASSERT_CALL(sleqp_aug_jacobian_set_iterate(jacobian, iterate));

  ASSERT_CALL(sleqp_dual_estimation_data_create(&estimation_data, problem));

  ASSERT_CALL(sleqp_dual_estimation_compute(estimation_data,
                                            iterate,
                                            NULL,
                                            jacobian));

  ASSERT_CALL(sleqp_sparse_vector_fprintf(iterate->cons_dual, stdout));

  ck_assert_int_eq(cons_dual->dim, 2);
  ck_assert_int_eq(cons_dual->nnz, 1);

  ck_assert_int_eq(cons_dual->dim, 2);

  ck_assert(sleqp_eq(cons_dual->data[0], 0.41421356237309515));

  ASSERT_CALL(sleqp_dual_estimation_data_free(&estimation_data));

  ASSERT_CALL(sleqp_aug_jacobian_free(&jacobian));

}
END_TEST

void constrained_teardown()
{
  ASSERT_CALL(sleqp_sparse_vector_free(&cauchy_direction));

  ASSERT_CALL(sleqp_cauchy_data_free(&cauchy_data));

  ASSERT_CALL(sleqp_lpi_free(&lp_interface));

  ASSERT_CALL(sleqp_iterate_free(&iterate));

  ASSERT_CALL(sleqp_problem_free(&problem));

  quadconsfunc_teardown();
}

Suite* constrained_test_suite()
{
  Suite *suite;
  TCase *tc_cons;

  suite = suite_create("Constrained tests");

  tc_cons = tcase_create("Constrained Cauchy direction");

  tcase_add_checked_fixture(tc_cons,
                            constrained_setup,
                            constrained_teardown);

  tcase_add_test(tc_cons, test_active_set);
  tcase_add_test(tc_cons, test_dual_variable);

  suite_add_tcase(suite, tc_cons);

  return suite;
}

int main()
{
  int num_fails;
  Suite* suite;
  SRunner* srunner;

  suite = constrained_test_suite();
  srunner = srunner_create(suite);

  srunner_set_fork_status(srunner, CK_NOFORK);
  srunner_run_all(srunner, CK_NORMAL);

  num_fails = srunner_ntests_failed(srunner);

  srunner_free(srunner);

  return (num_fails > 0) ? EXIT_FAILURE : EXIT_SUCCESS;
}
