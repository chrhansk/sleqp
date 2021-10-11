#include <stdlib.h>
#include <check.h>

#include "cmp.h"
#include "gauss_newton.h"
#include "problem.h"
#include "working_step.h"

#include "aug_jac/unconstrained_aug_jac.h"

#include "test_common.h"

#include "linear_lsq.h"

SleqpProblem* problem;

SleqpWorkingStep* working_step;

SleqpGaussNewtonSolver* gauss_newton_solver;

SleqpIterate* iterate;

SleqpAugJac* aug_jac;

SleqpSparseVec* direction;

void gauss_newton_setup()
{
  linear_lsq_setup();

  ASSERT_CALL(sleqp_problem_create_simple(&problem,
                                          linear_lsq_func,
                                          linear_lsq_params,
                                          linear_lsq_var_lb,
                                          linear_lsq_var_ub,
                                          linear_lsq_cons_lb,
                                          linear_lsq_cons_ub));

  ASSERT_CALL(sleqp_working_step_create(&working_step,
                                        problem,
                                        linear_lsq_params));

  ASSERT_CALL(sleqp_gauss_newton_solver_create(&gauss_newton_solver,
                                               problem,
                                               working_step,
                                               linear_lsq_params));

  ASSERT_CALL(sleqp_iterate_create(&iterate, problem, linear_lsq_initial));

  ASSERT_CALL(sleqp_unconstrained_aug_jac_create(&aug_jac,
                                                 problem));

  ASSERT_CALL(sleqp_sparse_vector_create_empty(&direction, linear_lsq_num_variables));
}

void gauss_newton_teardown()
{
  ASSERT_CALL(sleqp_sparse_vector_free(&direction));

  ASSERT_CALL(sleqp_aug_jac_release(&aug_jac));

  ASSERT_CALL(sleqp_iterate_release(&iterate));

  ASSERT_CALL(sleqp_gauss_newton_solver_release(&gauss_newton_solver));

  ASSERT_CALL(sleqp_working_step_release(&working_step));

  ASSERT_CALL(sleqp_problem_release(&problem));

  linear_lsq_teardown();
}

// LSQR direction should point to optimum
START_TEST(test_solve)
{
  const double penalty_parameter = 1.;
  const double trust_radius = 100.;

  const double eps = sleqp_params_get(linear_lsq_params,
                                      SLEQP_PARAM_EPS);

  ASSERT_CALL(sleqp_gauss_newton_solver_set_iterate(gauss_newton_solver,
                                                    iterate,
                                                    aug_jac,
                                                    trust_radius,
                                                    penalty_parameter));

  ASSERT_CALL(sleqp_gauss_newton_solver_compute_step(gauss_newton_solver,
                                                     direction));

  ck_assert(sleqp_sparse_vector_eq(direction, linear_lsq_optimal, eps));
}
END_TEST

// Insufficient trust radius, solution should be on the boundary
START_TEST(test_solve_small)
{
  const double penalty_parameter = 1.;
  const double trust_radius = 1.;

  const double eps = sleqp_params_get(linear_lsq_params,
                                      SLEQP_PARAM_EPS);

  ASSERT_CALL(sleqp_gauss_newton_solver_set_iterate(gauss_newton_solver,
                                                    iterate,
                                                    aug_jac,
                                                    trust_radius,
                                                    penalty_parameter));

  ASSERT_CALL(sleqp_gauss_newton_solver_compute_step(gauss_newton_solver,
                                                     direction));

  ck_assert(sleqp_is_eq(sleqp_sparse_vector_norm(direction),
                        trust_radius,
                        eps));
}
END_TEST

Suite* gauss_newton_test_suite()
{
  Suite *suite;
  TCase *tc_solve;

  suite = suite_create("Gauss-Newton tests");

  tc_solve = tcase_create("Solution test");

  tcase_add_checked_fixture(tc_solve,
                            gauss_newton_setup,
                            gauss_newton_teardown);

  tcase_add_test(tc_solve, test_solve);

  tcase_add_test(tc_solve, test_solve_small);

  suite_add_tcase(suite, tc_solve);

  return suite;
}

TEST_MAIN(gauss_newton_test_suite)
