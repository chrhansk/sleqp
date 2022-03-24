#include <check.h>
#include <stdlib.h>

#include "iterate.h"
#include "polish.h"
#include "problem.h"
#include "working_set.h"

#include "test_common.h"
#include "zero_func.h"

const int num_variables   = 3;
const int num_constraints = 3;

SleqpFunc* func;

SleqpParams* params;

SleqpVec* var_lb;
SleqpVec* var_ub;

SleqpVec* cons_lb;
SleqpVec* cons_ub;

SleqpProblem* problem;

SleqpVec* primal;

SleqpIterate* iterate;

SleqpPolishing* polishing;

void
polish_setup()
{
  ASSERT_CALL(zero_func_create(&func, num_variables, num_constraints));

  ASSERT_CALL(sleqp_params_create(&params));

  const double zero_eps = sleqp_params_value(params, SLEQP_PARAM_ZERO_EPS);

  ASSERT_CALL(sleqp_vec_create_empty(&var_lb, num_variables));
  ASSERT_CALL(sleqp_vec_create_full(&var_ub, num_variables));

  {
    double ub[] = {2., 2., 2.};

    ASSERT_CALL(sleqp_vec_from_raw(var_ub, ub, num_variables, zero_eps));
  }

  ASSERT_CALL(sleqp_vec_create_empty(&cons_lb, num_constraints));
  ASSERT_CALL(sleqp_vec_create_full(&cons_ub, num_constraints));

  {
    double ub[] = {2., 2., 2.};

    ASSERT_CALL(sleqp_vec_from_raw(cons_ub, ub, num_constraints, zero_eps));
  }

  ASSERT_CALL(sleqp_problem_create_simple(&problem,
                                          func,
                                          params,
                                          var_lb,
                                          var_ub,
                                          cons_lb,
                                          cons_ub));

  ASSERT_CALL(sleqp_vec_create_full(&primal, num_variables));

  double values[] = {0., 1., 2.};

  {
    ASSERT_CALL(sleqp_vec_from_raw(primal, values, num_variables, zero_eps));
  }

  ASSERT_CALL(sleqp_iterate_create(&iterate, problem, primal));

  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);

  ASSERT_CALL(
    sleqp_vec_reserve(sleqp_iterate_vars_dual(iterate), num_variables));

  ASSERT_CALL(
    sleqp_vec_reserve(sleqp_iterate_cons_dual(iterate), num_constraints));

  ASSERT_CALL(sleqp_working_set_reset(working_set));

  {
    ASSERT_CALL(sleqp_vec_from_raw(sleqp_iterate_cons_val(iterate),
                                   values,
                                   num_constraints,
                                   zero_eps));
  }

  ASSERT_CALL(sleqp_polishing_create(&polishing, problem, params));
}

void
polish_teardown()
{
  ASSERT_CALL(sleqp_polishing_release(&polishing));

  ASSERT_CALL(sleqp_iterate_release(&iterate));

  ASSERT_CALL(sleqp_vec_free(&primal));

  ASSERT_CALL(sleqp_problem_release(&problem));

  ASSERT_CALL(sleqp_vec_free(&cons_ub));
  ASSERT_CALL(sleqp_vec_free(&cons_lb));

  ASSERT_CALL(sleqp_vec_free(&var_ub));
  ASSERT_CALL(sleqp_vec_free(&var_lb));

  ASSERT_CALL(sleqp_params_release(&params));

  ASSERT_CALL(sleqp_func_release(&func));
}

START_TEST(test_polish_inactive_vars)
{
  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);

  // zero dual but active => kept
  ASSERT_CALL(sleqp_working_set_add_var(working_set, 0, SLEQP_ACTIVE_LOWER));

  // zero dual and inactive => removed
  ASSERT_CALL(sleqp_working_set_add_var(working_set, 1, SLEQP_ACTIVE_BOTH));

  ASSERT_CALL(
    sleqp_polishing_polish(polishing, iterate, SLEQP_POLISHING_INACTIVE));

  ck_assert_int_eq(sleqp_working_set_size(working_set), 1);

  ck_assert_int_eq(sleqp_working_set_var_state(working_set, 0),
                   SLEQP_ACTIVE_LOWER);

  ck_assert_int_eq(sleqp_working_set_var_state(working_set, 1), SLEQP_INACTIVE);
}
END_TEST

START_TEST(test_polish_zero_dual_vars)
{
  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);

  ASSERT_CALL(sleqp_working_set_add_var(working_set, 0, SLEQP_ACTIVE_LOWER));

  ASSERT_CALL(sleqp_working_set_add_var(working_set, 2, SLEQP_ACTIVE_UPPER));

  ASSERT_CALL(sleqp_vec_push(sleqp_iterate_vars_dual(iterate), 0, -1.));

  ASSERT_CALL(
    sleqp_polishing_polish(polishing, iterate, SLEQP_POLISHING_ZERO_DUAL));

  ck_assert_int_eq(sleqp_working_set_size(working_set), 1);

  ck_assert_int_eq(sleqp_working_set_var_state(working_set, 0),
                   SLEQP_ACTIVE_LOWER);

  ck_assert_int_eq(sleqp_working_set_var_state(working_set, 1), SLEQP_INACTIVE);

  ck_assert_int_eq(sleqp_working_set_var_state(working_set, 2), SLEQP_INACTIVE);
}
END_TEST

START_TEST(test_polish_inactive_cons)
{
  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);

  // zero dual but active => kept
  ASSERT_CALL(sleqp_working_set_add_cons(working_set, 0, SLEQP_ACTIVE_LOWER));

  // zero dual and inactive => removed
  ASSERT_CALL(sleqp_working_set_add_cons(working_set, 1, SLEQP_ACTIVE_BOTH));

  ASSERT_CALL(
    sleqp_polishing_polish(polishing, iterate, SLEQP_POLISHING_INACTIVE));

  ck_assert_int_eq(sleqp_working_set_size(working_set), 1);

  ck_assert_int_eq(sleqp_working_set_cons_state(working_set, 0),
                   SLEQP_ACTIVE_LOWER);

  ck_assert_int_eq(sleqp_working_set_cons_state(working_set, 1),
                   SLEQP_INACTIVE);
}
END_TEST

START_TEST(test_polish_zero_dual_cons)
{
  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);

  ASSERT_CALL(sleqp_working_set_add_cons(working_set, 0, SLEQP_ACTIVE_LOWER));

  ASSERT_CALL(sleqp_working_set_add_cons(working_set, 2, SLEQP_ACTIVE_UPPER));

  ASSERT_CALL(sleqp_vec_push(sleqp_iterate_cons_dual(iterate), 0, -1.));

  ASSERT_CALL(
    sleqp_polishing_polish(polishing, iterate, SLEQP_POLISHING_ZERO_DUAL));

  ck_assert_int_eq(sleqp_working_set_size(working_set), 1);

  ck_assert_int_eq(sleqp_working_set_cons_state(working_set, 0),
                   SLEQP_ACTIVE_LOWER);

  ck_assert_int_eq(sleqp_working_set_cons_state(working_set, 1),
                   SLEQP_INACTIVE);

  ck_assert_int_eq(sleqp_working_set_cons_state(working_set, 2),
                   SLEQP_INACTIVE);
}
END_TEST

Suite*
polish_test_suite()
{
  Suite* suite;
  TCase* tc_vars;
  TCase* tc_cons;

  suite = suite_create("Polishing tests");

  tc_vars = tcase_create("Variable polishing");

  tcase_add_checked_fixture(tc_vars, polish_setup, polish_teardown);

  tcase_add_test(tc_vars, test_polish_zero_dual_vars);

  tcase_add_test(tc_vars, test_polish_inactive_vars);

  suite_add_tcase(suite, tc_vars);

  tc_cons = tcase_create("Constraint polishing");

  tcase_add_checked_fixture(tc_cons, polish_setup, polish_teardown);

  tcase_add_test(tc_cons, test_polish_zero_dual_cons);

  tcase_add_test(tc_cons, test_polish_inactive_cons);

  suite_add_tcase(suite, tc_cons);

  return suite;
}

TEST_MAIN(polish_test_suite)
