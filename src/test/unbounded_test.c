#include <stdlib.h>
#include <check.h>

#include "cmp.h"
#include "mem.h"
#include "solver.h"

#include "test_common.h"

double v;

const int num_variables = 1;
const int num_constraints = 0;

SleqpFunc* unbounded_func;

SleqpSparseVec* unbounded_initial;

SleqpSparseVec* unbounded_var_lb;
SleqpSparseVec* unbounded_var_ub;
SleqpSparseVec* unbounded_cons_lb;
SleqpSparseVec* unbounded_cons_ub;

static SLEQP_RETCODE unbounded_set(SleqpFunc* func,
                                   SleqpSparseVec* x,
                                   SLEQP_VALUE_REASON reason,
                                   int* func_grad_nnz,
                                   int* cons_val_nnz,
                                   int* cons_jac_nnz,
                                   void* func_data)
{
  assert(x->dim == 1);
  SLEQP_CALL(sleqp_sparse_vector_to_raw(x, &v));

  (*func_grad_nnz) = 1;
  (*cons_val_nnz) = 1;
  (*cons_jac_nnz) = 1;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE unbounded_val(SleqpFunc* func,
                                   double* func_val,
                                   void* func_data)
{
  (*func_val) = v;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE unbounded_grad(SleqpFunc* func,
                                    SleqpSparseVec* func_grad,
                                    void* func_data)
{
  SLEQP_CALL(sleqp_sparse_vector_push(func_grad, 0, 1.));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE unbounded_hess_prod(SleqpFunc* func,
                                         const double* func_dual,
                                         const SleqpSparseVec* direction,
                                         const SleqpSparseVec* cons_duals,
                                         SleqpSparseVec* product,
                                         void* func_data)
{
  SLEQP_CALL(sleqp_sparse_vector_clear(product));

  return SLEQP_OKAY;
}

void unbounded_setup()
{
  const double inf = sleqp_infinity();

  SleqpFuncCallbacks callbacks = {
    .set_value = unbounded_set,
    .func_val  = unbounded_val,
    .func_grad = unbounded_grad,
    .cons_val  = NULL,
    .cons_jac  = NULL,
    .hess_prod = unbounded_hess_prod,
    .func_free = NULL
  };

  ASSERT_CALL(sleqp_func_create(&unbounded_func,
                                &callbacks,
                                num_variables,
                                num_constraints,
                                NULL));

  ASSERT_CALL(sleqp_sparse_vector_create_empty(&unbounded_initial, num_variables));

  ASSERT_CALL(sleqp_sparse_vector_create_full(&unbounded_var_lb, num_variables));
  ASSERT_CALL(sleqp_sparse_vector_push(unbounded_var_lb, 0, -inf));
  ASSERT_CALL(sleqp_sparse_vector_create_full(&unbounded_var_ub, num_variables));
  ASSERT_CALL(sleqp_sparse_vector_push(unbounded_var_ub, 0, inf));

  ASSERT_CALL(sleqp_sparse_vector_create_empty(&unbounded_cons_lb, num_constraints));
  ASSERT_CALL(sleqp_sparse_vector_create_empty(&unbounded_cons_ub, num_constraints));
}

void unbounded_teardown()
{
  ASSERT_CALL(sleqp_sparse_vector_free(&unbounded_cons_ub));
  ASSERT_CALL(sleqp_sparse_vector_free(&unbounded_cons_lb));

  ASSERT_CALL(sleqp_sparse_vector_free(&unbounded_var_ub));
  ASSERT_CALL(sleqp_sparse_vector_free(&unbounded_var_lb));

  ASSERT_CALL(sleqp_sparse_vector_free(&unbounded_initial));

  ASSERT_CALL(sleqp_func_release(&unbounded_func));
}

START_TEST(test_unbounded_solve)
{
  SleqpParams* params;
  SleqpOptions* options;
  SleqpProblem* problem;
  SleqpSolver* solver;

  ASSERT_CALL(sleqp_params_create(&params));

  ASSERT_CALL(sleqp_options_create(&options));

  ASSERT_CALL(sleqp_problem_create_simple(&problem,
                                          unbounded_func,
                                          params,
                                          unbounded_var_lb,
                                          unbounded_var_ub,
                                          unbounded_cons_lb,
                                          unbounded_cons_ub));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  params,
                                  options,
                                  unbounded_initial,
                                  NULL));

  // 100 iterations should be plenty...
  ASSERT_CALL(sleqp_solver_solve(solver, 100, -1));

  SleqpIterate* solution_iterate;

  ASSERT_CALL(sleqp_solver_get_solution(solver,
                                        &solution_iterate));

  ck_assert_int_eq(sleqp_solver_get_status(solver), SLEQP_STATUS_UNBOUNDED);

  ASSERT_CALL(sleqp_solver_release(&solver));

  ASSERT_CALL(sleqp_problem_release(&problem));

  ASSERT_CALL(sleqp_options_release(&options));

  ASSERT_CALL(sleqp_params_release(&params));
}
END_TEST


Suite* unconstrained_test_suite()
{
  Suite *suite;
  TCase *tc_unbounded;

  suite = suite_create("Unbounded test");

  tc_unbounded = tcase_create("Unbounded solution test");

  tcase_add_checked_fixture(tc_unbounded,
                            unbounded_setup,
                            unbounded_teardown);

  tcase_add_test(tc_unbounded, test_unbounded_solve);
  suite_add_tcase(suite, tc_unbounded);

  return suite;
}

TEST_MAIN(unconstrained_test_suite)
