#include <stdlib.h>
#include <check.h>

#include "sleqp_cmp.h"
#include "sleqp_mem.h"
#include "sleqp_solver.h"

#include "test_common.h"

SleqpSparseVec* var_lb;
SleqpSparseVec* var_ub;
SleqpSparseVec* cons_lb;
SleqpSparseVec* cons_ub;
SleqpSparseVec* x;

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

static SLEQP_RETCODE func_set(SleqpSparseVec* x,
                              int num_variables,
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

static SLEQP_RETCODE func_eval(int num_variables,
                               SleqpSparseVec* cons_indices,
                               double* func_val,
                               SleqpSparseVec* func_grad,
                               SleqpSparseVec* cons_val,
                               SleqpSparseMatrix* cons_jac,
                               void* func_data)
{
  FuncData* data = (FuncData*) func_data;
  double* x = data->values;

  if(func_val)
  {
    (*func_val) = x[0]*x[3]*(x[0] + x[1]+ x[2]) + x[2];
  }

  if(cons_val)
  {
    SLEQP_CALL(sleqp_sparse_vector_push(cons_val,
                                        0,
                                        x[0]*x[1]*x[2]*x[3]));

    SLEQP_CALL(sleqp_sparse_vector_push(cons_val,
                                        1,
                                        sq(x[0]) + sq(x[1]) + sq(x[2]) +  sq(x[3])));
  }

  if(func_grad)
  {
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
  }

  if(cons_jac)
  {
    SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac,
                                        0,
                                        0,
                                        x[1]*x[2]*x[3]));

    SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac,
                                        1,
                                        0,
                                        2*x[0]));

    SLEQP_CALL(sleqp_sparse_matrix_add_column(cons_jac, 1));

    SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac,
                                        0,
                                        1,
                                        x[0]*x[2]*x[3]));

    SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac,
                                        1,
                                        1,
                                        2*x[1]));

    SLEQP_CALL(sleqp_sparse_matrix_add_column(cons_jac, 2));

    SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac,
                                        0,
                                        2,
                                        x[0]*x[1]*x[3]));

    SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac,
                                        1,
                                        2,
                                        2*x[2]));

    SLEQP_CALL(sleqp_sparse_matrix_add_column(cons_jac, 3));

    SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac,
                                        0,
                                        3,
                                        x[0]*x[1]*x[2]));

    SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac,
                                        1,
                                        3,
                                        2*x[3]));

  }


  return SLEQP_OKAY;
}

static SLEQP_RETCODE func_hess_prod(int num_variables,
                                    double* func_dual,
                                    SleqpSparseVec* direction,
                                    SleqpSparseVec* cons_duals,
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

  ASSERT_CALL(sleqp_calloc(&func_data->values, num_variables));
  ASSERT_CALL(sleqp_calloc(&func_data->duals, num_constraints));
  ASSERT_CALL(sleqp_calloc(&func_data->direction, num_variables));

  ASSERT_CALL(sleqp_func_create(&func,
                                func_set,
                                func_eval,
                                func_hess_prod,
                                num_variables,
                                func_data));

  ASSERT_CALL(sleqp_params_create(&params));

  ASSERT_CALL(sleqp_problem_create(&problem,
                                   func,
                                   params,
                                   var_lb,
                                   var_ub,
                                   cons_lb,
                                   cons_ub));
}

void constrained_teardown()
{
  ASSERT_CALL(sleqp_problem_free(&problem));

  ASSERT_CALL(sleqp_params_free(&params));

  ASSERT_CALL(sleqp_func_free(&func));

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

START_TEST(test_constrained_solve)
{
  SleqpSparseVec* expected_solution;

  SleqpParams* params;
  SleqpSolver* solver;

  ASSERT_CALL(sleqp_sparse_vector_create(&expected_solution, 4, 4));

  /*
  (0) = 1.000000e+00
  (1) = 4.742999e+00
  (2) = 3.821151e+00
  (3) = 1.379408e+00
  */

  ASSERT_CALL(sleqp_sparse_vector_push(expected_solution, 0, 1.));
  ASSERT_CALL(sleqp_sparse_vector_push(expected_solution, 1, 4.742999));
  ASSERT_CALL(sleqp_sparse_vector_push(expected_solution, 2, 3.821151));
  ASSERT_CALL(sleqp_sparse_vector_push(expected_solution, 3, 1.379408));


  ASSERT_CALL(sleqp_params_create(&params));

  ASSERT_CALL(sleqp_solver_create(&solver,
                                  problem,
                                  params,
                                  x));

  // 100 iterations should be plenty...
  ASSERT_CALL(sleqp_solver_solve(solver, 100, -1));

  SleqpIterate* solution_iterate;

  ASSERT_CALL(sleqp_solver_get_solution(solver,
                                        &solution_iterate));

  ck_assert_int_eq(sleqp_solver_get_status(solver), SLEQP_OPTIMAL);

  SleqpSparseVec* actual_solution = solution_iterate->x;

  ck_assert(sleqp_sparse_vector_eq(actual_solution,
                                   expected_solution,
                                   1e-6));

  ASSERT_CALL(sleqp_solver_free(&solver));

  ASSERT_CALL(sleqp_params_free(&params));

  ASSERT_CALL(sleqp_sparse_vector_free(&expected_solution));
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

  tcase_add_test(tc_cons, test_constrained_solve);
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
