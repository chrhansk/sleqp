#include <stdlib.h>
#include <check.h>

#include "aug_jacobian.h"
#include "cauchy.h"
#include "cmp.h"
#include "dual_estimation.h"
#include "mem.h"
#include "util.h"

#include "sparse/sparse_factorization_umfpack.h"

#include "lp/lpi.h"

#include "test_common.h"

#include "quadcons_fixture.h"

SleqpParams* params;
SleqpOptions* options;
SleqpProblem* problem;
SleqpIterate* iterate;
SleqpLPi* lp_interface;
SleqpCauchy* cauchy_data;

SleqpSparseVec* cauchy_direction;

void constrained_setup()
{
  quadconsfunc_setup();

  ASSERT_CALL(sleqp_params_create(&params));

  ASSERT_CALL(sleqp_options_create(&options));

  ASSERT_CALL(sleqp_problem_create_simple(&problem,
                                          quadconsfunc,
                                          params,
                                          quadconsfunc_var_lb,
                                          quadconsfunc_var_ub,
                                          quadconsfunc_cons_lb,
                                          quadconsfunc_cons_ub));

  ASSERT_CALL(sleqp_iterate_create(&iterate,
                                   problem,
                                   quadconsfunc_x));

  const int num_variables = sleqp_problem_num_variables(problem);
  const int num_constraints = sleqp_problem_num_constraints(problem);

  int num_lp_variables = num_variables + 2*num_constraints;
  int num_lp_constraints = num_constraints;

  ASSERT_CALL(sleqp_lpi_create_default_interface(&lp_interface,
                                                 num_lp_variables,
                                                 num_lp_constraints,
                                                 params,
                                                 options));

  ASSERT_CALL(sleqp_set_and_evaluate(problem,
                                     iterate,
                                     SLEQP_VALUE_REASON_NONE,
                                     NULL));

  ASSERT_CALL(sleqp_cauchy_create(&cauchy_data,
                                  problem,
                                  params,
                                  options,
                                  lp_interface));

  ASSERT_CALL(sleqp_sparse_vector_create(&cauchy_direction, 2, 2));

}

START_TEST(test_working_set)
{
  SleqpWorkingSet* working_set = sleqp_iterate_get_working_set(iterate);

  double trust_radius = 0.1, penalty_parameter = 1.;

  ASSERT_CALL(sleqp_cauchy_set_iterate(cauchy_data,
                                       iterate,
                                       trust_radius));

  ASSERT_CALL(sleqp_cauchy_solve(cauchy_data,
                                 sleqp_iterate_get_func_grad(iterate),
                                 penalty_parameter,
                                 SLEQP_CAUCHY_OBJECTIVE_TYPE_DEFAULT));

  ASSERT_CALL(sleqp_cauchy_get_working_set(cauchy_data,
                                           iterate));

  ck_assert_int_eq(sleqp_working_set_get_constraint_state(working_set, 0), SLEQP_INACTIVE);
  ck_assert_int_eq(sleqp_working_set_get_constraint_state(working_set, 1), SLEQP_ACTIVE_UPPER);
}
END_TEST

START_TEST(test_dual_variable)
{
  SleqpSparseVec* cons_dual = sleqp_iterate_get_cons_dual(iterate);

  SleqpSparseFactorization* factorization;
  SleqpAugJacobian* jacobian;
  SleqpDualEstimation* estimation_data;

  double trust_radius = 0.1, penalty_parameter = 1.;

  ASSERT_CALL(sleqp_cauchy_set_iterate(cauchy_data,
                                       iterate,
                                       trust_radius));

  ASSERT_CALL(sleqp_cauchy_solve(cauchy_data,
                                 sleqp_iterate_get_func_grad(iterate),
                                 penalty_parameter,
                                 SLEQP_CAUCHY_OBJECTIVE_TYPE_DEFAULT));

  ASSERT_CALL(sleqp_cauchy_get_working_set(cauchy_data,
                                           iterate));

  ASSERT_CALL(sleqp_cauchy_get_direction(cauchy_data, cauchy_direction));

  ASSERT_CALL(sleqp_sparse_factorization_create_default(&factorization,
                                                        params));

  ASSERT_CALL(sleqp_aug_jacobian_create(&jacobian,
                                        problem,
                                        params,
                                        factorization));

  ASSERT_CALL(sleqp_aug_jacobian_set_iterate(jacobian, iterate));

  ASSERT_CALL(sleqp_dual_estimation_create(&estimation_data, problem));

  ASSERT_CALL(sleqp_dual_estimation_compute(estimation_data,
                                            iterate,
                                            NULL,
                                            jacobian));

  ck_assert_int_eq(cons_dual->dim, 2);
  ck_assert_int_eq(cons_dual->nnz, 1);

  ck_assert_int_eq(cons_dual->dim, 2);

  ck_assert(sleqp_is_eq(cons_dual->data[0], 0.4142135623, 1e-8));

  ASSERT_CALL(sleqp_dual_estimation_free(&estimation_data));

  ASSERT_CALL(sleqp_aug_jacobian_release(&jacobian));

  ASSERT_CALL(sleqp_sparse_factorization_release(&factorization));

}
END_TEST

void constrained_teardown()
{
  ASSERT_CALL(sleqp_sparse_vector_free(&cauchy_direction));

  ASSERT_CALL(sleqp_cauchy_release(&cauchy_data));

  ASSERT_CALL(sleqp_lpi_free(&lp_interface));

  ASSERT_CALL(sleqp_iterate_release(&iterate));

  ASSERT_CALL(sleqp_problem_release(&problem));

  ASSERT_CALL(sleqp_options_release(&options));

  ASSERT_CALL(sleqp_params_release(&params));

  quadconsfunc_teardown();
}

Suite* constrained_cauchy_test_suite()
{
  Suite *suite;
  TCase *tc_cons;

  suite = suite_create("Constrained tests");

  tc_cons = tcase_create("Constrained Cauchy direction");

  tcase_add_checked_fixture(tc_cons,
                            constrained_setup,
                            constrained_teardown);

  tcase_add_test(tc_cons, test_working_set);
  tcase_add_test(tc_cons, test_dual_variable);

  suite_add_tcase(suite, tc_cons);

  return suite;
}

TEST_MAIN(constrained_cauchy_test_suite)
