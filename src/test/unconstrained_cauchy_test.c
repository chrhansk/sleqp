
#include <stdlib.h>
#include <check.h>

#include "cmp.h"
#include "params.h"

#include "cauchy/unconstrained_cauchy.h"

#include "test_common.h"
#include "zero_func.h"

SleqpParams* params;

const int num_variables = 2;
const int num_constraints = 0;

const double objective = 25.;

SleqpFunc* func;

SleqpSparseVec* var_lb;
SleqpSparseVec* var_ub;
SleqpSparseVec* cons_lb;
SleqpSparseVec* cons_ub;
SleqpSparseVec* primal;
SleqpSparseVec* grad;

SleqpProblem* problem;
SleqpIterate* iterate;

SleqpSparseVec* direction;

SleqpCauchy* cauchy;

void unconstrained_setup()
{
  const double inf = sleqp_infinity();

  double neg_inf_vals[] = {inf, inf};
  double pos_inf_vals[] = {inf, inf};

  ASSERT_CALL(sleqp_params_create(&params));

  ASSERT_CALL(zero_func_create(&func,
                               num_variables,
                               num_constraints));

  ASSERT_CALL(sleqp_sparse_vector_create_full(&var_lb,
                                               num_variables));

  ASSERT_CALL(sleqp_sparse_vector_create_full(&var_ub,
                                              num_variables));

  const double zero_eps = sleqp_params_get(params,
                                           SLEQP_PARAM_ZERO_EPS);

  ASSERT_CALL(sleqp_sparse_vector_from_raw(var_lb,
                                           neg_inf_vals,
                                           num_variables,
                                           zero_eps));

  ASSERT_CALL(sleqp_sparse_vector_from_raw(var_ub,
                                           pos_inf_vals,
                                           num_variables,
                                           zero_eps));

  ASSERT_CALL(sleqp_sparse_vector_create_empty(&cons_lb,
                                               num_constraints));

  ASSERT_CALL(sleqp_sparse_vector_create_empty(&cons_ub,
                                               num_constraints));

  ASSERT_CALL(sleqp_sparse_vector_create_full(&primal,
                                              num_variables));

  double primal_vals[] = {1., 1.};

  ASSERT_CALL(sleqp_sparse_vector_from_raw(primal,
                                           primal_vals,
                                           num_variables,
                                           zero_eps));

  ASSERT_CALL(sleqp_sparse_vector_create_full(&grad,
                                              num_variables));

  double grad_vals[] = {1., -1.};

  ASSERT_CALL(sleqp_sparse_vector_from_raw(grad,
                                           grad_vals,
                                           num_variables,
                                           zero_eps));

  ASSERT_CALL(sleqp_problem_create_simple(&problem,
                                          func,
                                          params,
                                          var_lb,
                                          var_ub,
                                          cons_lb,
                                          cons_ub));

  ASSERT_CALL(sleqp_iterate_create(&iterate, problem, primal));

  ASSERT_CALL(sleqp_sparse_vector_copy(grad,
                                       sleqp_iterate_get_func_grad(iterate)));

  ASSERT_CALL(sleqp_iterate_set_func_val(iterate, objective));

  ASSERT_CALL(sleqp_sparse_vector_create_empty(&direction,
                                               num_variables));

  ASSERT_CALL(sleqp_unconstrained_cauchy_create(&cauchy,
                                                problem,
                                                params));
}

void unconstrained_teardown()
{
  ASSERT_CALL(sleqp_cauchy_release(&cauchy));

  ASSERT_CALL(sleqp_sparse_vector_free(&direction));

  ASSERT_CALL(sleqp_iterate_release(&iterate));

  ASSERT_CALL(sleqp_problem_release(&problem));

  ASSERT_CALL(sleqp_sparse_vector_free(&grad));
  ASSERT_CALL(sleqp_sparse_vector_free(&primal));

  ASSERT_CALL(sleqp_sparse_vector_free(&cons_ub));
  ASSERT_CALL(sleqp_sparse_vector_free(&cons_lb));

  ASSERT_CALL(sleqp_sparse_vector_free(&var_ub));
  ASSERT_CALL(sleqp_sparse_vector_free(&var_lb));

  ASSERT_CALL(sleqp_func_release(&func));

  ASSERT_CALL(sleqp_params_release(&params));
}

START_TEST(test_solve)
{
  const double trust_radius = 100.;

  ASSERT_CALL(sleqp_cauchy_set_iterate(cauchy,
                                       iterate,
                                       trust_radius));

  ASSERT_CALL(sleqp_cauchy_solve(cauchy,
                                 grad,
                                 1.,
                                 SLEQP_CAUCHY_OBJECTIVE_TYPE_DEFAULT));

  ASSERT_CALL(sleqp_cauchy_get_working_set(cauchy,
                                           iterate));

  SleqpWorkingSet* working_set = sleqp_iterate_get_working_set(iterate);

  ck_assert_int_eq(sleqp_working_set_get_variable_state(working_set, 0),
                   SLEQP_INACTIVE);

  ck_assert_int_eq(sleqp_working_set_get_variable_state(working_set, 1),
                   SLEQP_INACTIVE);

  ASSERT_CALL(sleqp_cauchy_get_direction(cauchy, direction));

  ck_assert_int_eq(sleqp_sparse_vector_value_at(direction, 0),
                   -trust_radius);

  ck_assert_int_eq(sleqp_sparse_vector_value_at(direction, 1),
                   trust_radius);

  const double eps = sleqp_params_get(params, SLEQP_PARAM_EPS);

  double actual_objective;

  ASSERT_CALL(sleqp_cauchy_get_objective_value(cauchy,
                                               &actual_objective));

  double inner_product;

  ASSERT_CALL(sleqp_sparse_vector_dot(direction,
                                      grad,
                                      &inner_product));

  const double expected_objective = sleqp_iterate_get_func_val(iterate) + inner_product;

  ck_assert(sleqp_is_eq(actual_objective,
                        expected_objective,
                        eps));

  ASSERT_CALL(sleqp_cauchy_get_dual_estimation(cauchy, iterate));

  SleqpSparseVec* vars_dual = sleqp_iterate_get_vars_dual(iterate);

  ck_assert(vars_dual->nnz == 0);
}
END_TEST

Suite* unconstrained_cauchy_test_suite()
{
  Suite *suite;
  TCase *tc_cons;

  suite = suite_create("Unconstrained Cauchy suite");

  tc_cons = tcase_create("Unconstrained Cauchy tests");

  tcase_add_checked_fixture(tc_cons,
                            unconstrained_setup,
                            unconstrained_teardown);

  tcase_add_test(tc_cons, test_solve);

  suite_add_tcase(suite, tc_cons);

  return suite;
}

TEST_MAIN(unconstrained_cauchy_test_suite)
