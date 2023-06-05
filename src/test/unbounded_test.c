#include <check.h>
#include <stdlib.h>

#include "cmp.h"
#include "mem.h"
#include "solver.h"

#include "test_common.h"

double v;

const int num_variables   = 1;
const int num_constraints = 0;

SleqpFunc* unbounded_func;

SleqpVec* unbounded_initial;

SleqpVec* unbounded_var_lb;
SleqpVec* unbounded_var_ub;
SleqpVec* unbounded_cons_lb;
SleqpVec* unbounded_cons_ub;

static SLEQP_RETCODE
unbounded_set(SleqpFunc* func,
              SleqpVec* x,
              SLEQP_VALUE_REASON reason,
              bool* reject,
              void* func_data)
{
  assert(x->dim == 1);
  SLEQP_CALL(sleqp_vec_to_raw(x, &v));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
unbounded_obj_val(SleqpFunc* func, double* obj_val, void* func_data)
{
  (*obj_val) = v;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
unbounded_obj_grad(SleqpFunc* func, SleqpVec* obj_grad, void* func_data)
{
  SLEQP_CALL(sleqp_vec_push(obj_grad, 0, 1.));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
unbounded_hess_prod(SleqpFunc* func,
                    const SleqpVec* direction,
                    const SleqpVec* cons_duals,
                    SleqpVec* product,
                    void* func_data)
{
  SLEQP_CALL(sleqp_vec_clear(product));

  return SLEQP_OKAY;
}

void
unbounded_setup()
{
  const double inf = sleqp_infinity();

  SleqpFuncCallbacks callbacks = {.set_value = unbounded_set,
                                  .obj_val   = unbounded_obj_val,
                                  .obj_grad  = unbounded_obj_grad,
                                  .cons_val  = NULL,
                                  .cons_jac  = NULL,
                                  .hess_prod = unbounded_hess_prod,
                                  .func_free = NULL};

  ASSERT_CALL(sleqp_func_create(&unbounded_func,
                                &callbacks,
                                num_variables,
                                num_constraints,
                                NULL));

  ASSERT_CALL(sleqp_vec_create_empty(&unbounded_initial, num_variables));

  ASSERT_CALL(sleqp_vec_create_full(&unbounded_var_lb, num_variables));
  ASSERT_CALL(sleqp_vec_push(unbounded_var_lb, 0, -inf));
  ASSERT_CALL(sleqp_vec_create_full(&unbounded_var_ub, num_variables));
  ASSERT_CALL(sleqp_vec_push(unbounded_var_ub, 0, inf));

  ASSERT_CALL(sleqp_vec_create_empty(&unbounded_cons_lb, num_constraints));
  ASSERT_CALL(sleqp_vec_create_empty(&unbounded_cons_ub, num_constraints));
}

void
unbounded_teardown()
{
  ASSERT_CALL(sleqp_vec_free(&unbounded_cons_ub));
  ASSERT_CALL(sleqp_vec_free(&unbounded_cons_lb));

  ASSERT_CALL(sleqp_vec_free(&unbounded_var_ub));
  ASSERT_CALL(sleqp_vec_free(&unbounded_var_lb));

  ASSERT_CALL(sleqp_vec_free(&unbounded_initial));

  ASSERT_CALL(sleqp_func_release(&unbounded_func));
}

START_TEST(test_unbounded_solve)
{
  SleqpSettings* settings;
  SleqpProblem* problem;
  SleqpSolver* solver;

  ASSERT_CALL(sleqp_settings_create(&settings));

  ASSERT_CALL(sleqp_problem_create_simple(&problem,
                                          unbounded_func,
                                          unbounded_var_lb,
                                          unbounded_var_ub,
                                          unbounded_cons_lb,
                                          unbounded_cons_ub,
                                          settings));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  unbounded_initial,
                                  NULL));

  // 100 iterations should be plenty...
  ASSERT_CALL(sleqp_solver_solve(solver, 100, -1));

  SleqpIterate* solution_iterate;

  ASSERT_CALL(sleqp_solver_solution(solver, &solution_iterate));

  ck_assert_int_eq(sleqp_solver_status(solver), SLEQP_STATUS_UNBOUNDED);

  ASSERT_CALL(sleqp_solver_release(&solver));

  ASSERT_CALL(sleqp_problem_release(&problem));

  ASSERT_CALL(sleqp_settings_release(&settings));
}
END_TEST

Suite*
unconstrained_test_suite()
{
  Suite* suite;
  TCase* tc_unbounded;

  suite = suite_create("Unbounded test");

  tc_unbounded = tcase_create("Unbounded solution test");

  tcase_add_checked_fixture(tc_unbounded, unbounded_setup, unbounded_teardown);

  tcase_add_test(tc_unbounded, test_unbounded_solve);
  suite_add_tcase(suite, tc_unbounded);

  return suite;
}

TEST_MAIN(unconstrained_test_suite)
