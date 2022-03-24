#include <check.h>
#include <stdlib.h>
#include <unistd.h>

#include "test_common.h"
#include "zero_func.h"

#include "options.h"
#include "problem.h"
#include "solver.h"

const int num_variables   = 1;
const int num_constraints = 1;

const double delay_time = 1e-4;

double value;

void
delay()
{
  clock_t start = clock();

  // spin until time limit reached
  while (true)
  {
    for (volatile int i = 0; i < 10000; ++i)
      ;

    clock_t current = clock();

    double elapsed = ((double)(current - start)) / CLOCKS_PER_SEC;

    if (elapsed >= delay_time)
    {
      break;
    }
  }
}

static SLEQP_RETCODE
delay_func_set(SleqpFunc* func,
               SleqpVec* x,
               SLEQP_VALUE_REASON reason,
               bool* reject,
               int* obj_grad_nnz,
               int* cons_val_nnz,
               int* cons_jac_nnz,
               void* func_data)
{
  *obj_grad_nnz = 1;
  *cons_val_nnz = 0;
  *cons_jac_nnz = 0;

  sleqp_vec_value_at(x, value);

  delay();

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
delay_func_obj_val(SleqpFunc* func, double* obj_val, void* func_data)
{
  *obj_val = value;

  delay();

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
delay_func_obj_grad(SleqpFunc* func, SleqpVec* obj_grad, void* func_data)
{
  SLEQP_CALL(sleqp_vec_push(obj_grad, 0, 1.));

  delay();

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
delay_func_cons_val(SleqpFunc* func, SleqpVec* cons_val, void* func_data)
{
  return SLEQP_OKAY;
}

static SLEQP_RETCODE
delay_func_cons_jac(SleqpFunc* func,
                    SleqpSparseMatrix* cons_jac,
                    void* func_data)
{
  return SLEQP_OKAY;
}

static SLEQP_RETCODE
delay_func_hess_prod(SleqpFunc* func,
                     const double* obj_dual,
                     const SleqpVec* direction,
                     const SleqpVec* cons_duals,
                     SleqpVec* result,
                     void* func_data)
{
  return SLEQP_OKAY;
}

SleqpFunc* func;

SleqpParams* params;
SleqpOptions* options;

SleqpVec* var_lb;
SleqpVec* var_ub;

SleqpVec* cons_lb;
SleqpVec* cons_ub;

SleqpProblem* problem;

SleqpVec* primal;

void
time_limit_setup()
{
  SleqpFuncCallbacks callbacks = {.set_value = delay_func_set,
                                  .obj_val   = delay_func_obj_val,
                                  .obj_grad  = delay_func_obj_grad,
                                  .cons_val  = delay_func_cons_val,
                                  .cons_jac  = delay_func_cons_jac,
                                  .hess_prod = delay_func_hess_prod,
                                  .func_free = NULL};

  ASSERT_CALL(
    sleqp_func_create(&func, &callbacks, num_variables, num_constraints, NULL));

  ASSERT_CALL(sleqp_params_create(&params));

  ASSERT_CALL(sleqp_options_create(&options));

  ASSERT_CALL(sleqp_vec_create_empty(&var_lb, num_variables));
  ASSERT_CALL(sleqp_vec_create_empty(&var_ub, num_variables));

  ASSERT_CALL(sleqp_vec_create_empty(&cons_lb, num_constraints));
  ASSERT_CALL(sleqp_vec_create_empty(&cons_ub, num_constraints));

  ASSERT_CALL(sleqp_problem_create_simple(&problem,
                                          func,
                                          params,
                                          var_lb,
                                          var_ub,
                                          cons_lb,
                                          cons_ub));

  ASSERT_CALL(sleqp_vec_create_empty(&primal, num_variables));
}

void
time_limit_teardown()
{
  ASSERT_CALL(sleqp_vec_free(&primal));

  ASSERT_CALL(sleqp_problem_release(&problem));

  ASSERT_CALL(sleqp_vec_free(&cons_ub));
  ASSERT_CALL(sleqp_vec_free(&cons_lb));

  ASSERT_CALL(sleqp_vec_free(&var_ub));
  ASSERT_CALL(sleqp_vec_free(&var_lb));

  ASSERT_CALL(sleqp_options_release(&options));

  ASSERT_CALL(sleqp_params_release(&params));

  ASSERT_CALL(sleqp_func_release(&func));
}

const double time_limit = 1e-6;

START_TEST(test_solve)
{
  SleqpSolver* solver;

  ASSERT_CALL(
    sleqp_solver_create(&solver, problem, params, options, primal, NULL));

  ASSERT_CALL(sleqp_solver_solve(solver, SLEQP_NONE, time_limit));

  ck_assert_int_eq(sleqp_solver_status(solver), SLEQP_STATUS_ABORT_TIME);

  ASSERT_CALL(sleqp_solver_release(&solver));
}
END_TEST

Suite*
time_limit_test_suite()
{
  Suite* suite;
  TCase* tc_solve;

  suite = suite_create("Time limit tests");

  tc_solve = tcase_create("Solution test");

  tcase_add_checked_fixture(tc_solve, time_limit_setup, time_limit_teardown);

  tcase_add_test(tc_solve, test_solve);

  suite_add_tcase(suite, tc_solve);

  return suite;
}

TEST_MAIN(time_limit_test_suite)
