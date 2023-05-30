#include <check.h>
#include <stdlib.h>

#include "cmp.h"
#include "params.h"

#include "cauchy/unconstrained_cauchy.h"

#include "test_common.h"
#include "zero_func.h"

SleqpParams* params;

const int num_variables   = 2;
const int num_constraints = 0;

const double objective = 25.;

SleqpFunc* func;

SleqpVec* var_lb;
SleqpVec* var_ub;
SleqpVec* cons_lb;
SleqpVec* cons_ub;
SleqpVec* primal;
SleqpVec* grad;

SleqpProblem* problem;
SleqpIterate* iterate;

SleqpVec* direction;

SleqpCauchy* cauchy;

void
unconstrained_setup()
{
  const double inf = sleqp_infinity();

  ASSERT_CALL(sleqp_params_create(&params));

  ASSERT_CALL(zero_func_create(&func, num_variables, num_constraints));

  ASSERT_CALL(sleqp_vec_create_full(&var_lb, num_variables));
  ASSERT_CALL(sleqp_vec_fill(var_lb, -inf));

  ASSERT_CALL(sleqp_vec_create_full(&var_ub, num_variables));
  ASSERT_CALL(sleqp_vec_fill(var_ub, -inf));

  const double zero_eps = sleqp_params_value(params, SLEQP_PARAM_ZERO_EPS);

  ASSERT_CALL(sleqp_vec_create_empty(&cons_lb, num_constraints));

  ASSERT_CALL(sleqp_vec_create_empty(&cons_ub, num_constraints));

  ASSERT_CALL(sleqp_vec_create_full(&primal, num_variables));

  double primal_vals[] = {1., 1.};

  ASSERT_CALL(
    sleqp_vec_set_from_raw(primal, primal_vals, num_variables, zero_eps));

  ASSERT_CALL(sleqp_vec_create_full(&grad, num_variables));

  double grad_vals[] = {1., -1.};

  ASSERT_CALL(sleqp_vec_set_from_raw(grad, grad_vals, num_variables, zero_eps));

  ASSERT_CALL(sleqp_problem_create_simple(&problem,
                                          func,
                                          params,
                                          var_lb,
                                          var_ub,
                                          cons_lb,
                                          cons_ub));

  ASSERT_CALL(sleqp_iterate_create(&iterate, problem, primal));

  ASSERT_CALL(sleqp_vec_copy(grad, sleqp_iterate_obj_grad(iterate)));

  ASSERT_CALL(sleqp_iterate_set_obj_val(iterate, objective));

  ASSERT_CALL(sleqp_vec_create_empty(&direction, num_variables));

  ASSERT_CALL(sleqp_unconstrained_cauchy_create(&cauchy, problem, params));
}

void
unconstrained_teardown()
{
  ASSERT_CALL(sleqp_cauchy_release(&cauchy));

  ASSERT_CALL(sleqp_vec_free(&direction));

  ASSERT_CALL(sleqp_iterate_release(&iterate));

  ASSERT_CALL(sleqp_problem_release(&problem));

  ASSERT_CALL(sleqp_vec_free(&grad));
  ASSERT_CALL(sleqp_vec_free(&primal));

  ASSERT_CALL(sleqp_vec_free(&cons_ub));
  ASSERT_CALL(sleqp_vec_free(&cons_lb));

  ASSERT_CALL(sleqp_vec_free(&var_ub));
  ASSERT_CALL(sleqp_vec_free(&var_lb));

  ASSERT_CALL(sleqp_func_release(&func));

  ASSERT_CALL(sleqp_params_release(&params));
}

START_TEST(test_solve)
{
  const double trust_radius = 100.;

  ASSERT_CALL(sleqp_cauchy_set_iterate(cauchy, iterate, trust_radius));

  ASSERT_CALL(
    sleqp_cauchy_solve(cauchy, grad, 1., SLEQP_CAUCHY_OBJTYPE_DEFAULT));

  ASSERT_CALL(sleqp_cauchy_working_set(cauchy, iterate));

  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);

  ck_assert_int_eq(sleqp_working_set_var_state(working_set, 0), SLEQP_INACTIVE);

  ck_assert_int_eq(sleqp_working_set_var_state(working_set, 1), SLEQP_INACTIVE);

  ASSERT_CALL(sleqp_cauchy_lp_step(cauchy, direction));

  ck_assert_int_eq(sleqp_vec_value_at(direction, 0), -trust_radius);

  ck_assert_int_eq(sleqp_vec_value_at(direction, 1), trust_radius);

  const double eps = sleqp_params_value(params, SLEQP_PARAM_EPS);

  double actual_objective;

  ASSERT_CALL(sleqp_cauchy_obj_val(cauchy, &actual_objective));

  double inner_product;

  ASSERT_CALL(sleqp_vec_dot(direction, grad, &inner_product));

  const double expected_objective
    = sleqp_iterate_obj_val(iterate) + inner_product;

  ck_assert(sleqp_is_eq(actual_objective, expected_objective, eps));

  ASSERT_CALL(sleqp_cauchy_estimate_duals(cauchy,
                                          sleqp_iterate_working_set(iterate),
                                          sleqp_iterate_cons_dual(iterate),
                                          sleqp_iterate_vars_dual(iterate)));

  SleqpVec* vars_dual = sleqp_iterate_vars_dual(iterate);

  ck_assert(vars_dual->nnz == 0);
}
END_TEST

Suite*
unconstrained_cauchy_test_suite()
{
  Suite* suite;
  TCase* tc_cons;

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
