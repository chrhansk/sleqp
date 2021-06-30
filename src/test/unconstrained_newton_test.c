#include <stdlib.h>
#include <check.h>

#include "sleqp.h"
#include "sleqp_aug_jacobian.h"
#include "sleqp_cauchy.h"
#include "sleqp_cmp.h"
#include "sleqp_dual_estimation.h"
#include "sleqp_mem.h"
#include "sleqp_newton.h"

#include "sparse/sleqp_sparse_factorization_umfpack.h"

#include "test_common.h"

#include "quadfunc_fixture.h"

SleqpParams* params;
SleqpOptions* options;
SleqpProblem* problem;
SleqpIterate* iterate;
SleqpLPi* lp_interface;
SleqpCauchy* cauchy_data;

const double tolerance = 1e-8;

void newton_setup()
{
  quadfunc_setup();

  ASSERT_CALL(sleqp_params_create(&params));

  ASSERT_CALL(sleqp_options_create(&options));

  ASSERT_CALL(sleqp_problem_create_simple(&problem,
                                          quadfunc,
                                          params,
                                          quadfunc_var_lb,
                                          quadfunc_var_ub,
                                          quadfunc_cons_lb,
                                          quadfunc_cons_ub));

  const int num_variables = sleqp_problem_num_variables(problem);
  const int num_constraints = sleqp_problem_num_constraints(problem);

  ASSERT_CALL(sleqp_iterate_create(&iterate,
                                   problem,
                                   quadfunc_x));

  SleqpWorkingSet* working_set = sleqp_iterate_get_working_set(iterate);
  ASSERT_CALL(sleqp_working_set_reset(working_set));

  int num_lp_variables = num_variables + 2*num_constraints;
  int num_lp_constraints = num_constraints;

  ASSERT_CALL(sleqp_lpi_create_default_interface(&lp_interface,
                                                 num_lp_variables,
                                                 num_lp_constraints,
                                                 params));

  ASSERT_CALL(sleqp_set_and_evaluate(problem, iterate, SLEQP_VALUE_REASON_INIT));

  ASSERT_CALL(sleqp_cauchy_create(&cauchy_data,
                                  problem,
                                  params,
                                  options,
                                  lp_interface));
}

void newton_teardown()
{
  ASSERT_CALL(sleqp_cauchy_free(&cauchy_data));

  ASSERT_CALL(sleqp_lpi_free(&lp_interface));

  ASSERT_CALL(sleqp_iterate_release(&iterate));

  ASSERT_CALL(sleqp_problem_release(&problem));

  ASSERT_CALL(sleqp_options_release(&options));

  ASSERT_CALL(sleqp_params_release(&params));

  quadfunc_teardown();
}

START_TEST(newton_wide_step)
{
  SleqpNewtonData* newton_data;

  SleqpSparseVec* expected_step;
  SleqpSparseVec* actual_step;

  SleqpSparseFactorization* factorization;
  SleqpAugJacobian* jacobian;

  const int num_variables = sleqp_problem_num_variables(problem);

  ASSERT_CALL(sleqp_sparse_vector_create(&expected_step, num_variables, 2));

  ASSERT_CALL(sleqp_sparse_vector_push(expected_step, 0, -1.));
  ASSERT_CALL(sleqp_sparse_vector_push(expected_step, 1, -2.));

  ASSERT_CALL(sleqp_sparse_vector_create(&actual_step, num_variables, 0));

  double penalty_parameter = 1.;
  double trust_radius = 10.;

  ASSERT_CALL(sleqp_sparse_factorization_create_default(&factorization,
                                                        params));

  // create with empty active set
  ASSERT_CALL(sleqp_aug_jacobian_create(&jacobian,
                                        problem,
                                        params,
                                        factorization));

  ASSERT_CALL(sleqp_aug_jacobian_set_iterate(jacobian, iterate));

  ASSERT_CALL(sleqp_newton_data_create(&newton_data,
                                       problem,
                                       params,
                                       options));

  ASSERT_CALL(sleqp_newton_set_iterate(newton_data,
                                       iterate,
                                       jacobian,
                                       trust_radius,
                                       penalty_parameter));

  // we use the default (empty) active set for the Newton step,
  // trust region size should be large to ensure that
  // the solution is that of the unrestricted step
  ASSERT_CALL(sleqp_newton_compute_step(newton_data,
                                        sleqp_iterate_get_cons_dual(iterate),
                                        actual_step));

  ck_assert(sleqp_sparse_vector_eq(expected_step, actual_step, tolerance));

  ASSERT_CALL(sleqp_newton_data_release(&newton_data));

  ASSERT_CALL(sleqp_aug_jacobian_release(&jacobian));

  ASSERT_CALL(sleqp_sparse_factorization_release(&factorization));

  ASSERT_CALL(sleqp_sparse_vector_free(&actual_step));

  ASSERT_CALL(sleqp_sparse_vector_free(&expected_step));
}
END_TEST

START_TEST(newton_small_step)
{
  SleqpNewtonData* newton_data;

  SleqpSparseVec* expected_step;
  SleqpSparseVec* actual_step;

  SleqpSparseFactorization* factorization;
  SleqpAugJacobian* jacobian;

  const int num_variables = sleqp_problem_num_variables(problem);

  ASSERT_CALL(sleqp_sparse_vector_create(&expected_step, num_variables, 2));

  ASSERT_CALL(sleqp_sparse_vector_push(expected_step, 0, -0.44721359549995793));
  ASSERT_CALL(sleqp_sparse_vector_push(expected_step, 1, -0.89442719099991586));

  ASSERT_CALL(sleqp_sparse_vector_create(&actual_step, num_variables, 0));

  double penalty_parameter = 1.;
  double trust_radius = 1.;

  ASSERT_CALL(sleqp_sparse_factorization_create_default(&factorization,
                                                        params));

  // create with empty active set
  ASSERT_CALL(sleqp_aug_jacobian_create(&jacobian,
                                        problem,
                                        params,
                                        factorization));

  ASSERT_CALL(sleqp_aug_jacobian_set_iterate(jacobian, iterate));

  ASSERT_CALL(sleqp_newton_data_create(&newton_data,
                                       problem,
                                       params,
                                       options));

  ASSERT_CALL(sleqp_newton_set_iterate(newton_data,
                                       iterate,
                                       jacobian,
                                       trust_radius,
                                       penalty_parameter));

  // we use the default (empty) active set for the Newton step,
  // trust region size should be so small that
  // the solution is on the boundary of the feasible set
  ASSERT_CALL(sleqp_newton_compute_step(newton_data,
                                        sleqp_iterate_get_cons_dual(iterate),
                                        actual_step));

  ck_assert(sleqp_sparse_vector_eq(expected_step, actual_step, tolerance));

  ASSERT_CALL(sleqp_newton_data_release(&newton_data));

  ASSERT_CALL(sleqp_aug_jacobian_release(&jacobian));

  ASSERT_CALL(sleqp_sparse_factorization_release(&factorization));

  ASSERT_CALL(sleqp_sparse_vector_free(&actual_step));

  ASSERT_CALL(sleqp_sparse_vector_free(&expected_step));
}
END_TEST

Suite* newton_test_suite()
{
  Suite *suite;
  TCase *tc_cons;

  suite = suite_create("Unconstrained newton step tests");

  tc_cons = tcase_create("Newton step");

  tcase_add_checked_fixture(tc_cons,
                            newton_setup,
                            newton_teardown);

  tcase_add_test(tc_cons, newton_wide_step);

  tcase_add_test(tc_cons, newton_small_step);

  suite_add_tcase(suite, tc_cons);

  return suite;
}

TEST_MAIN(newton_test_suite)
