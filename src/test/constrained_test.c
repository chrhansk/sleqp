#include <stdlib.h>
#include <check.h>

#include "sleqp_cmp.h"
#include "sleqp_mem.h"
#include "sleqp_solver.h"
#include "sleqp_util.h"

#include "test_common.h"

SleqpSparseVec* var_lb;
SleqpSparseVec* var_ub;
SleqpSparseVec* cons_lb;
SleqpSparseVec* cons_ub;
SleqpSparseVec* x;

SleqpParams* params;
SleqpOptions* options;

SleqpSparseVec* expected_solution;

const int num_variables = 4;
const int num_constraints = 2;

typedef struct FuncData
{
  double* values;
  double* duals;
  double* direction;

} FuncData;

static double sq(double x)
{
  return x*x;
}

static SLEQP_RETCODE func_set(SleqpFunc* func,
                              SleqpSparseVec* x,
                              SLEQP_VALUE_REASON reason,
                              int* func_grad_nnz,
                              int* cons_val_nnz,
                              int* cons_jac_nnz,
                              void* func_data)
{
  FuncData* data = (FuncData*) func_data;

  SLEQP_CALL(sleqp_sparse_vector_to_raw(x, data->values));

  *func_grad_nnz = num_variables;
  *cons_val_nnz = num_constraints;
  *cons_jac_nnz = num_constraints * num_variables;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE func_val(SleqpFunc* func,
                              double* func_val,
                              void* func_data)
{
  FuncData* data = (FuncData*) func_data;
  double* x = data->values;

  (*func_val) = x[0]*x[3]*(x[0] + x[1]+ x[2]) + x[2];

  return SLEQP_OKAY;
}

static SLEQP_RETCODE func_grad(SleqpFunc* func,
                               SleqpSparseVec* func_grad,
                               void* func_data)
{
  FuncData* data = (FuncData*) func_data;
  double* x = data->values;

  SLEQP_CALL(sleqp_sparse_vector_push(func_grad,
                                      0,
                                      (x[0] + x[1] + x[2])*x[3] + x[0]*x[3]));

  SLEQP_CALL(sleqp_sparse_vector_push(func_grad,
                                      1,
                                      x[0]*x[3]));

  SLEQP_CALL(sleqp_sparse_vector_push(func_grad,
                                      2,
                                      x[0]*x[3] + 1));

  SLEQP_CALL(sleqp_sparse_vector_push(func_grad,
                                      3,
                                      (x[0] + x[1] + x[2])*x[0]));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE cons_val(SleqpFunc* func,
                              const SleqpSparseVec* cons_indices,
                              SleqpSparseVec* cons_val,
                              void* func_data)
{
  FuncData* data = (FuncData*) func_data;
  double* x = data->values;

  SLEQP_CALL(sleqp_sparse_vector_push(cons_val,
                                      0,
                                      x[0]*x[1]*x[2]*x[3]));

  SLEQP_CALL(sleqp_sparse_vector_push(cons_val,
                                      1,
                                      sq(x[0]) + sq(x[1]) + sq(x[2]) +  sq(x[3])));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE cons_jac(SleqpFunc* func,
                              const SleqpSparseVec* cons_indices,
                              SleqpSparseMatrix* cons_jac,
                              void* func_data)
{
  FuncData* data = (FuncData*) func_data;
  double* x = data->values;

  SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac,
                                      0,
                                      0,
                                      x[1]*x[2]*x[3]));

  SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac,
                                      1,
                                      0,
                                      2*x[0]));

  SLEQP_CALL(sleqp_sparse_matrix_push_column(cons_jac, 1));

  SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac,
                                      0,
                                      1,
                                      x[0]*x[2]*x[3]));

  SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac,
                                      1,
                                      1,
                                      2*x[1]));

  SLEQP_CALL(sleqp_sparse_matrix_push_column(cons_jac, 2));

  SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac,
                                      0,
                                      2,
                                      x[0]*x[1]*x[3]));

  SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac,
                                      1,
                                      2,
                                      2*x[2]));

  SLEQP_CALL(sleqp_sparse_matrix_push_column(cons_jac, 3));

  SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac,
                                      0,
                                      3,
                                      x[0]*x[1]*x[2]));

  SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac,
                                      1,
                                      3,
                                      2*x[3]));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE func_hess_prod(SleqpFunc* func,
                                    const double* func_dual,
                                    const SleqpSparseVec* direction,
                                    const SleqpSparseVec* cons_duals,
                                    SleqpSparseVec* product,
                                    void* func_data)
{
  FuncData* data = (FuncData*) func_data;

  SLEQP_CALL(sleqp_sparse_vector_to_raw(cons_duals, data->duals));

  SLEQP_CALL(sleqp_sparse_vector_to_raw(direction, data->direction));

  double* x = data->values;

  double* duals = data->duals;
  double* dir = data->direction;

  double f_dual = func_dual ? *func_dual : 0.;

  SLEQP_CALL(sleqp_sparse_vector_reserve(product, num_variables));

  {
    double v = 0.;

    v += (2*x[3]*f_dual + 2*duals[1])*dir[0];
    v += (x[3]*f_dual + x[2]*x[3]*duals[0])*dir[1];
    v += (x[3]*f_dual + x[1]*x[3]*duals[0])*dir[2];
    v += ((2*x[0] +x[1] + x[2])*f_dual + x[1]*x[2]*duals[0])*dir[3];

    SLEQP_CALL(sleqp_sparse_vector_push(product, 0, v));
  }

  {
    double v = 0.;

    v += (x[3]*f_dual + x[2]*x[3]*duals[0])*dir[0];
    v += (2*duals[1])*dir[1];
    v += (x[0]*x[3]*duals[0])*dir[2];
    v += (x[0]*f_dual + (x[0]*x[2])*duals[0])*dir[3];

    SLEQP_CALL(sleqp_sparse_vector_push(product, 1, v));
  }

  {
    double v = 0.;

    v += (x[3]*f_dual + x[1]*x[3]*duals[0])*dir[0];
    v += (x[0]*x[3]*duals[0])*dir[1];
    v += (2*duals[1])*dir[2];
    v += (x[0]*f_dual + x[0]*x[1]*duals[0])*dir[3];

    SLEQP_CALL(sleqp_sparse_vector_push(product, 2, v));
  }

  {
    double v = 0.;

    v += ((2*x[0] + x[1] + x[2])*f_dual + x[1]*x[2]*duals[0])*dir[0];
    v += (x[0]*f_dual + x[0]*x[2]*duals[0])*dir[1];
    v += (x[0]*f_dual + x[0]*x[1]*duals[0])*dir[2];
    v += (2*duals[1])*dir[3];

    SLEQP_CALL(sleqp_sparse_vector_push(product, 3, v));
  }

  return SLEQP_OKAY;
}

FuncData* func_data;
SleqpFunc* func;

SleqpParams* params;
SleqpProblem* problem;

void constrained_setup()
{
  const double inf = sleqp_infinity();

  ASSERT_CALL(sleqp_sparse_vector_create(&var_lb,
                                         num_variables,
                                         num_variables));

  ASSERT_CALL(sleqp_sparse_vector_create(&var_ub,
                                         num_variables,
                                         num_variables));

  for(int i = 0; i < num_variables; ++i)
  {
    ASSERT_CALL(sleqp_sparse_vector_push(var_lb, i, 1.));
    ASSERT_CALL(sleqp_sparse_vector_push(var_ub, i, 5.));
  }

  ASSERT_CALL(sleqp_sparse_vector_create(&cons_lb,
                                         num_constraints,
                                         num_constraints));

  ASSERT_CALL(sleqp_sparse_vector_create(&cons_ub,
                                         num_constraints,
                                         num_constraints));

  ASSERT_CALL(sleqp_sparse_vector_push(cons_lb, 0, 25.));
  ASSERT_CALL(sleqp_sparse_vector_push(cons_lb, 1, 40.));

  ASSERT_CALL(sleqp_sparse_vector_push(cons_ub, 0, inf));
  ASSERT_CALL(sleqp_sparse_vector_push(cons_ub, 1, 40.));

  ASSERT_CALL(sleqp_sparse_vector_create(&x,
                                         num_variables,
                                         num_variables));

  ASSERT_CALL(sleqp_sparse_vector_push(x, 0, 1.));
  ASSERT_CALL(sleqp_sparse_vector_push(x, 1, 5.));
  ASSERT_CALL(sleqp_sparse_vector_push(x, 2, 5.));
  ASSERT_CALL(sleqp_sparse_vector_push(x, 3, 1.));

  ASSERT_CALL(sleqp_malloc(&func_data));

  ASSERT_CALL(sleqp_alloc_array(&func_data->values, num_variables));
  ASSERT_CALL(sleqp_alloc_array(&func_data->duals, num_constraints));
  ASSERT_CALL(sleqp_alloc_array(&func_data->direction, num_variables));

  SleqpFuncCallbacks callbacks = {
    .set_value = func_set,
    .func_val  = func_val,
    .func_grad = func_grad,
    .cons_val  = cons_val,
    .cons_jac  = cons_jac,
    .hess_prod = func_hess_prod,
    .func_free = NULL
  };

  ASSERT_CALL(sleqp_func_create(&func,
                                &callbacks,
                                num_variables,
                                num_constraints,
                                func_data));

  ASSERT_CALL(sleqp_params_create(&params));

  ASSERT_CALL(sleqp_options_create(&options));

  ASSERT_CALL(sleqp_problem_create(&problem,
                                   func,
                                   params,
                                   var_lb,
                                   var_ub,
                                   cons_lb,
                                   cons_ub));

  ASSERT_CALL(sleqp_sparse_vector_create(&expected_solution, 4, 4));

  ASSERT_CALL(sleqp_sparse_vector_push(expected_solution, 0, 1.));
  ASSERT_CALL(sleqp_sparse_vector_push(expected_solution, 1, 4.742999));
  ASSERT_CALL(sleqp_sparse_vector_push(expected_solution, 2, 3.821151));
  ASSERT_CALL(sleqp_sparse_vector_push(expected_solution, 3, 1.379408));
}

void constrained_teardown()
{
  ASSERT_CALL(sleqp_sparse_vector_free(&expected_solution));

  ASSERT_CALL(sleqp_problem_free(&problem));

  ASSERT_CALL(sleqp_options_release(&options));

  ASSERT_CALL(sleqp_params_release(&params));

  ASSERT_CALL(sleqp_func_release(&func));

  sleqp_free(&func_data->direction);
  sleqp_free(&func_data->duals);
  sleqp_free(&func_data->values);
  sleqp_free(&func_data);

  ASSERT_CALL(sleqp_sparse_vector_free(&x));
  ASSERT_CALL(sleqp_sparse_vector_free(&cons_ub));
  ASSERT_CALL(sleqp_sparse_vector_free(&cons_lb));
  ASSERT_CALL(sleqp_sparse_vector_free(&var_ub));
  ASSERT_CALL(sleqp_sparse_vector_free(&var_lb));
}

START_TEST(test_solve)
{
  SleqpSolver* solver;

  ASSERT_CALL(sleqp_options_set_deriv_check(options,
                                            SLEQP_DERIV_CHECK_FIRST | SLEQP_DERIV_CHECK_SECOND_EXHAUSTIVE));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  params,
                                  options,
                                  x,
                                  NULL));

  // 100 iterations should be plenty...
  ASSERT_CALL(sleqp_solver_solve(solver, 100, -1));

  SleqpIterate* solution_iterate;

  ASSERT_CALL(sleqp_solver_get_solution(solver,
                                        &solution_iterate));

  ck_assert_int_eq(sleqp_solver_get_status(solver), SLEQP_OPTIMAL);

  SleqpSparseVec* actual_solution = sleqp_iterate_get_primal(solution_iterate);

  ck_assert(sleqp_sparse_vector_eq(actual_solution,
                                   expected_solution,
                                   1e-6));

  ASSERT_CALL(sleqp_solver_release(&solver));
}
END_TEST

START_TEST(test_sr1_solve)
{
  ASSERT_CALL(sleqp_options_set_deriv_check(options,
                                            SLEQP_DERIV_CHECK_FIRST));

  ASSERT_CALL(sleqp_options_set_hessian_eval(options,
                                             SLEQP_HESSIAN_EVAL_SR1));

  SleqpSolver* solver;

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  params,
                                  options,
                                  x,
                                  NULL));

  // 100 iterations should be plenty...
  ASSERT_CALL(sleqp_solver_solve(solver, 100, -1));

  SleqpIterate* solution_iterate;

  ASSERT_CALL(sleqp_solver_get_solution(solver,
                                        &solution_iterate));

  ck_assert_int_eq(sleqp_solver_get_status(solver), SLEQP_OPTIMAL);

  SleqpSparseVec* actual_solution = sleqp_iterate_get_primal(solution_iterate);

  ck_assert(sleqp_sparse_vector_eq(actual_solution,
                                   expected_solution,
                                   1e-6));

  ASSERT_CALL(sleqp_solver_release(&solver));
}
END_TEST

START_TEST(test_bfgs_solve_no_sizing)
{
  ASSERT_CALL(sleqp_options_set_deriv_check(options,
                                            SLEQP_DERIV_CHECK_FIRST));

  ASSERT_CALL(sleqp_options_set_bfgs_sizing(options,
                                            SLEQP_BFGS_SIZING_NONE));

  ASSERT_CALL(sleqp_options_set_hessian_eval(options,
                                             SLEQP_HESSIAN_EVAL_DAMPED_BFGS));

  SleqpSolver* solver;

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  params,
                                  options,
                                  x,
                                  NULL));

  // 100 iterations should be plenty...
  ASSERT_CALL(sleqp_solver_solve(solver, 100, -1));

  SleqpIterate* solution_iterate;

  ASSERT_CALL(sleqp_solver_get_solution(solver,
                                        &solution_iterate));

  ck_assert_int_eq(sleqp_solver_get_status(solver), SLEQP_OPTIMAL);

  SleqpSparseVec* actual_solution = sleqp_iterate_get_primal(solution_iterate);

  ck_assert(sleqp_sparse_vector_eq(actual_solution,
                                   expected_solution,
                                   1e-6));

  ASSERT_CALL(sleqp_solver_release(&solver));
}
END_TEST

START_TEST(test_bfgs_solve_centered_ol_sizing)
{
  ASSERT_CALL(sleqp_options_set_deriv_check(options,
                                            SLEQP_DERIV_CHECK_FIRST));

  ASSERT_CALL(sleqp_options_set_bfgs_sizing(options,
                                            SLEQP_BFGS_SIZING_CENTERED_OL));

  ASSERT_CALL(sleqp_options_set_hessian_eval(options,
                                             SLEQP_HESSIAN_EVAL_DAMPED_BFGS));

  SleqpSolver* solver;

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  params,
                                  options,
                                  x,
                                  NULL));

  // 100 iterations should be plenty...
  ASSERT_CALL(sleqp_solver_solve(solver, 100, -1));

  SleqpIterate* solution_iterate;

  ASSERT_CALL(sleqp_solver_get_solution(solver,
                                        &solution_iterate));

  ck_assert_int_eq(sleqp_solver_get_status(solver), SLEQP_OPTIMAL);

  SleqpSparseVec* actual_solution = sleqp_iterate_get_primal(solution_iterate);

  ck_assert(sleqp_sparse_vector_eq(actual_solution,
                                   expected_solution,
                                   1e-6));

  ASSERT_CALL(sleqp_solver_release(&solver));
}
END_TEST

START_TEST(test_unscaled_solve)
{
  ASSERT_CALL(sleqp_options_set_deriv_check(options,
                                            SLEQP_DERIV_CHECK_FIRST));

  SleqpSolver* solver;

  SleqpScalingData* scaling_data;

  ASSERT_CALL(sleqp_scaling_create(&scaling_data,
                                   problem->num_variables,
                                   problem->num_constraints));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  params,
                                  options,
                                  x,
                                  scaling_data));

  // 100 iterations should be plenty...
  ASSERT_CALL(sleqp_solver_solve(solver, 100, -1));

  SleqpIterate* solution_iterate;

  ASSERT_CALL(sleqp_solver_get_solution(solver,
                                        &solution_iterate));

  ck_assert_int_eq(sleqp_solver_get_status(solver), SLEQP_OPTIMAL);

  SleqpSparseVec* actual_solution = sleqp_iterate_get_primal(solution_iterate);

  ck_assert(sleqp_sparse_vector_eq(actual_solution,
                                   expected_solution,
                                   1e-6));

  ASSERT_CALL(sleqp_solver_release(&solver));

  ASSERT_CALL(sleqp_scaling_release(&scaling_data));
}
END_TEST

START_TEST(test_scaled_solve)
{
  ASSERT_CALL(sleqp_options_set_deriv_check(options,
                                            SLEQP_DERIV_CHECK_FIRST | SLEQP_DERIV_CHECK_SECOND_EXHAUSTIVE));

  SleqpSolver* solver;

  SleqpScalingData* scaling_data;

  ASSERT_CALL(sleqp_scaling_create(&scaling_data,
                                   num_variables,
                                   num_constraints));

  ASSERT_CALL(sleqp_scaling_set_func_weight(scaling_data, 2));

  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling_data, 0, -5));
  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling_data, 1, 5));

  ASSERT_CALL(sleqp_scaling_set_cons_weight(scaling_data, 0, -1));
  ASSERT_CALL(sleqp_scaling_set_cons_weight(scaling_data, 1, 2));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  params,
                                  options,
                                  x,
                                  scaling_data));

  // 100 iterations should be plenty...
  ASSERT_CALL(sleqp_solver_solve(solver, 1000, -1));

  SleqpIterate* solution_iterate;

  ASSERT_CALL(sleqp_solver_get_solution(solver,
                                        &solution_iterate));

  ck_assert_int_eq(sleqp_solver_get_status(solver), SLEQP_OPTIMAL);

  SleqpSparseVec* actual_solution = sleqp_iterate_get_primal(solution_iterate);

  ck_assert(sleqp_sparse_vector_eq(actual_solution,
                                   expected_solution,
                                   1e-6));

  ASSERT_CALL(sleqp_solver_release(&solver));

  ASSERT_CALL(sleqp_scaling_release(&scaling_data));
}
END_TEST

START_TEST(test_scaled_sr1_solve)
{
  ASSERT_CALL(sleqp_options_set_deriv_check(options,
                                            SLEQP_DERIV_CHECK_FIRST));

  SleqpSolver* solver;

  SleqpScalingData* scaling_data;

  ASSERT_CALL(sleqp_options_set_hessian_eval(options,
                                             SLEQP_HESSIAN_EVAL_SR1));

  ASSERT_CALL(sleqp_scaling_create(&scaling_data,
                                   problem->num_variables,
                                   problem->num_constraints));

  ASSERT_CALL(sleqp_scaling_set_func_weight(scaling_data, 2));

  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling_data, 0, -5));
  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling_data, 1, 5));

  ASSERT_CALL(sleqp_scaling_set_cons_weight(scaling_data, 0, -1));
  ASSERT_CALL(sleqp_scaling_set_cons_weight(scaling_data, 1, 2));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  params,
                                  options,
                                  x,
                                  scaling_data));

  // 100 iterations should be plenty...
  ASSERT_CALL(sleqp_solver_solve(solver, 1000, -1));

  SleqpIterate* solution_iterate;

  ASSERT_CALL(sleqp_solver_get_solution(solver,
                                        &solution_iterate));

  ck_assert_int_eq(sleqp_solver_get_status(solver), SLEQP_OPTIMAL);

  SleqpSparseVec* actual_solution = sleqp_iterate_get_primal(solution_iterate);

  ck_assert(sleqp_sparse_vector_eq(actual_solution,
                                   expected_solution,
                                   1e-6));

  ASSERT_CALL(sleqp_solver_release(&solver));

  ASSERT_CALL(sleqp_scaling_release(&scaling_data));
}
END_TEST

START_TEST(test_scaled_bfgs_solve)
{
  ASSERT_CALL(sleqp_options_set_deriv_check(options,
                                            SLEQP_DERIV_CHECK_FIRST));

  SleqpSolver* solver;

  SleqpScalingData* scaling_data;

  ASSERT_CALL(sleqp_options_set_hessian_eval(options,
                                             SLEQP_HESSIAN_EVAL_DAMPED_BFGS));

  ASSERT_CALL(sleqp_scaling_create(&scaling_data,
                                   problem->num_variables,
                                   problem->num_constraints));

  ASSERT_CALL(sleqp_scaling_set_func_weight(scaling_data, 2));

  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling_data, 0, -5));
  ASSERT_CALL(sleqp_scaling_set_var_weight(scaling_data, 1, 5));

  ASSERT_CALL(sleqp_scaling_set_cons_weight(scaling_data, 0, -1));
  ASSERT_CALL(sleqp_scaling_set_cons_weight(scaling_data, 1, 2));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  params,
                                  options,
                                  x,
                                  scaling_data));

  // 100 iterations should be plenty...
  ASSERT_CALL(sleqp_solver_solve(solver, 1000, -1));

  SleqpIterate* solution_iterate;

  ASSERT_CALL(sleqp_solver_get_solution(solver,
                                        &solution_iterate));

  ck_assert_int_eq(sleqp_solver_get_status(solver), SLEQP_OPTIMAL);

  SleqpSparseVec* actual_solution = sleqp_iterate_get_primal(solution_iterate);

  ck_assert(sleqp_sparse_vector_eq(actual_solution,
                                   expected_solution,
                                   1e-6));

  ASSERT_CALL(sleqp_solver_release(&solver));

  ASSERT_CALL(sleqp_scaling_release(&scaling_data));
}
END_TEST

START_TEST(test_auto_scaled_solve)
{
  ASSERT_CALL(sleqp_options_set_deriv_check(options,
                                            SLEQP_DERIV_CHECK_FIRST));

  SleqpSolver* solver;

  SleqpScalingData* scaling_data;

  SleqpIterate* iterate;

  const double eps = sleqp_params_get_eps(params);

  ASSERT_CALL(sleqp_iterate_create(&iterate,
                                   problem,
                                   x));

  ASSERT_CALL(sleqp_set_and_evaluate(problem, iterate, SLEQP_VALUE_REASON_NONE));

  ASSERT_CALL(sleqp_scaling_create(&scaling_data,
                                   problem->num_variables,
                                   problem->num_constraints));

  ASSERT_CALL(sleqp_func_scaling_from_gradient(scaling_data,
                                               sleqp_iterate_get_func_grad(iterate),
                                               eps));

  ASSERT_CALL(sleqp_scaling_from_cons_jac(scaling_data,
                                          sleqp_iterate_get_cons_jac(iterate),
                                          eps));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  params,
                                  options,
                                  x,
                                  scaling_data));

  // 100 iterations should be plenty...
  ASSERT_CALL(sleqp_solver_solve(solver, 100, -1));

  SleqpIterate* solution_iterate;

  ASSERT_CALL(sleqp_solver_get_solution(solver,
                                        &solution_iterate));

  ck_assert_int_eq(sleqp_solver_get_status(solver), SLEQP_OPTIMAL);

  SleqpSparseVec* actual_solution = sleqp_iterate_get_primal(solution_iterate);

  ck_assert(sleqp_sparse_vector_eq(actual_solution,
                                   expected_solution,
                                   1e-6));

  ASSERT_CALL(sleqp_solver_release(&solver));

  ASSERT_CALL(sleqp_scaling_release(&scaling_data));

  ASSERT_CALL(sleqp_iterate_release(&iterate));
}
END_TEST

Suite* constrained_test_suite()
{
  Suite *suite;
  TCase *tc_cons;

  suite = suite_create("Constrained tests");

  tc_cons = tcase_create("Constrained solution test");

  tcase_add_checked_fixture(tc_cons,
                            constrained_setup,
                            constrained_teardown);

  tcase_add_test(tc_cons, test_solve);

  tcase_add_test(tc_cons, test_sr1_solve);

  tcase_add_test(tc_cons, test_bfgs_solve_no_sizing);

  tcase_add_test(tc_cons, test_bfgs_solve_centered_ol_sizing);

  tcase_add_test(tc_cons, test_unscaled_solve);

  tcase_add_test(tc_cons, test_scaled_solve);

  tcase_add_test(tc_cons, test_scaled_sr1_solve);

  // tcase_add_test(tc_cons, test_scaled_bfgs_solve);

  tcase_add_test(tc_cons, test_auto_scaled_solve);

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
