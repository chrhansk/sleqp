#include <stdlib.h>
#include <check.h>

#include "test_common.h"

#include "sleqp.h"
#include "sleqp_cmp.h"
#include "sleqp_func.h"
#include "sleqp_iterate.h"
#include "sleqp_mem.h"
#include "sleqp_newton.h"
#include "sleqp_problem.h"

#include "sparse/sleqp_sparse_factorization_umfpack.h"

static const int num_variables = 2;
static const int num_constraints = 1;

SleqpFunc* linquadfunc;

SleqpSparseVec* linquadfunc_var_lb;
SleqpSparseVec* linquadfunc_var_ub;
SleqpSparseVec* linquadfunc_cons_lb;
SleqpSparseVec* linquadfunc_cons_ub;
SleqpSparseVec* linquadfunc_x;

SleqpParams* params;
SleqpProblem* problem;
SleqpIterate* iterate;

typedef struct LinQuadFuncData
{
  double* x;
} LinQuadFuncData;

LinQuadFuncData* func_data;

static inline double square(double v)
{
  return v*v;
}

SLEQP_RETCODE linquadfunc_set(SleqpFunc* func,
                              SleqpSparseVec* x,
                              SLEQP_VALUE_REASON reason,
                              int* func_grad_nnz,
                              int* cons_val_nnz,
                              int* cons_jac_nnz,
                              void* func_data)
{
  *func_grad_nnz = 2;
  *cons_val_nnz = 1;
  *cons_jac_nnz = 1;

  LinQuadFuncData* data = (LinQuadFuncData*) func_data;

  data->x[0] = 0;
  data->x[1] = 0;

  int k_x = 0;

  while(k_x < x->nnz)
  {
    data->x[x->indices[k_x]] = x->data[k_x];

    ++k_x;
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE linquadfunc_val(SleqpFunc* func,
                              double* func_val,
                              void* func_data)
{
  LinQuadFuncData* data = (LinQuadFuncData*) func_data;

  *func_val = square(data->x[0]) + square(data->x[1]);

  return SLEQP_OKAY;
}

SLEQP_RETCODE linquadfunc_grad(SleqpFunc* func,
                               SleqpSparseVec* func_grad,
                               void* func_data)
{
  LinQuadFuncData* data = (LinQuadFuncData*) func_data;

  assert(func_grad->dim == 2);
  assert(func_grad->nnz_max >= 2);

  func_grad->nnz = 0;

  SLEQP_CALL(sleqp_sparse_vector_push(func_grad,
                                      0,
                                      2.*data->x[0]));

  SLEQP_CALL(sleqp_sparse_vector_push(func_grad,
                                      1,
                                      2.*data->x[1]));

  return SLEQP_OKAY;
}

SLEQP_RETCODE linquadfunc_cons_val(SleqpFunc* func,
                                   const SleqpSparseVec* cons_indices,
                                   SleqpSparseVec* cons_val,
                                   void* func_data)
{
  LinQuadFuncData* data = (LinQuadFuncData*) func_data;

  SLEQP_CALL(sleqp_sparse_vector_push(cons_val,
                                      0,
                                      data->x[1]));

  return SLEQP_OKAY;
}

SLEQP_RETCODE linquadfunc_cons_jac(SleqpFunc* func,
                                   const SleqpSparseVec* cons_indices,
                                   SleqpSparseMatrix* cons_jac,
                                   void* func_data)
{
  assert(sleqp_sparse_matrix_get_nnz(cons_jac) == 0);
  assert(sleqp_sparse_matrix_get_nnz_max(cons_jac) >= 1);

  SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac, 0, 1, 1.));

  assert(sleqp_sparse_matrix_valid(cons_jac));

  return SLEQP_OKAY;
}

SLEQP_RETCODE linquadfunc_hess_prod(SleqpFunc* func,
                                    const double* func_dual,
                                    const SleqpSparseVec* direction,
                                    const SleqpSparseVec* cons_duals,
                                    SleqpSparseVec* result,
                                    void* func_data)
{
  if(func_dual)
  {
    double total_value = 2.* (*func_dual);

    SLEQP_CALL(sleqp_sparse_vector_copy(direction, result));

    SLEQP_CALL(sleqp_sparse_vector_scale(result, total_value));
  }

  return SLEQP_OKAY;
}

void newton_setup()
{
  // setup function

  ASSERT_CALL(sleqp_malloc(&func_data));

  ASSERT_CALL(sleqp_alloc_array(&func_data->x, 2));

  SleqpFuncCallbacks callbacks = {
    .set_value = linquadfunc_set,
    .func_val  = linquadfunc_val,
    .func_grad = linquadfunc_grad,
    .cons_val  = linquadfunc_cons_val,
    .cons_jac  = linquadfunc_cons_jac,
    .hess_prod = linquadfunc_hess_prod,
    .func_free = NULL
  };

  ASSERT_CALL(sleqp_func_create(&linquadfunc,
                                &callbacks,
                                num_variables,
                                num_constraints,
                                func_data));

  double inf = sleqp_infinity();

  ASSERT_CALL(sleqp_sparse_vector_create_full(&linquadfunc_var_lb,
                                              num_variables));

  ASSERT_CALL(sleqp_sparse_vector_push(linquadfunc_var_lb, 0, -inf));
  ASSERT_CALL(sleqp_sparse_vector_push(linquadfunc_var_lb, 1, -inf));

  ASSERT_CALL(sleqp_sparse_vector_create_full(&linquadfunc_var_ub,
                                              num_variables));

  ASSERT_CALL(sleqp_sparse_vector_push(linquadfunc_var_ub, 0, inf));
  ASSERT_CALL(sleqp_sparse_vector_push(linquadfunc_var_ub, 1, inf));

  ASSERT_CALL(sleqp_sparse_vector_create_full(&linquadfunc_cons_lb,
                                              num_constraints));

  ASSERT_CALL(sleqp_sparse_vector_push(linquadfunc_cons_lb, 0, 2.));

  ASSERT_CALL(sleqp_sparse_vector_create_full(&linquadfunc_cons_ub,
                                              num_constraints));

  ASSERT_CALL(sleqp_sparse_vector_push(linquadfunc_cons_ub, 0, inf));

  ASSERT_CALL(sleqp_sparse_vector_create_full(&linquadfunc_x,
                                              num_variables));


  ASSERT_CALL(sleqp_sparse_vector_push(linquadfunc_x, 0, 1.));

  ASSERT_CALL(sleqp_sparse_vector_push(linquadfunc_x, 1, 2.));

  ASSERT_CALL(sleqp_params_create(&params));

  ASSERT_CALL(sleqp_problem_create(&problem,
                                   linquadfunc,
                                   params,
                                   linquadfunc_var_lb,
                                   linquadfunc_var_ub,
                                   linquadfunc_cons_lb,
                                   linquadfunc_cons_ub));

  ASSERT_CALL(sleqp_iterate_create(&iterate,
                                   problem,
                                   linquadfunc_x));

  ASSERT_CALL(sleqp_set_and_evaluate(problem, iterate, SLEQP_VALUE_REASON_NONE));

  SleqpWorkingSet* working_set = sleqp_iterate_get_working_set(iterate);
  ASSERT_CALL(sleqp_working_set_reset(working_set));

  // set the cons state manually...

  ASSERT_CALL(sleqp_working_set_add_constraint(working_set,
                                               0,
                                               SLEQP_ACTIVE_LOWER));
}

START_TEST(newton_constrained_step)
{
  SleqpSparseVec* actual_step;
  SleqpSparseVec* expected_step;

  SleqpParams* params;
  SleqpOptions* options;

  SleqpNewtonData* newton_data;

  SleqpSparseFactorization* factorization;
  SleqpAugJacobian* jacobian;

  double penalty_parameter = 1.;
  double trust_radius = 10.;

  int num_variables = problem->num_variables;

  ASSERT_CALL(sleqp_sparse_vector_create(&actual_step, num_variables, 0));
  ASSERT_CALL(sleqp_sparse_vector_create(&expected_step, num_variables, 1));

  ASSERT_CALL(sleqp_sparse_vector_push(expected_step, 0, -1.));

  ASSERT_CALL(sleqp_params_create(&params));

  ASSERT_CALL(sleqp_options_create(&options));

  ASSERT_CALL(sleqp_sparse_factorization_create_default(&factorization,
                                                        params));

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

  ASSERT_CALL(sleqp_newton_compute_step(newton_data,
                                        sleqp_iterate_get_cons_dual(iterate),
                                        actual_step));

  const double tolerance = 1e-8;

  ck_assert(sleqp_sparse_vector_eq(actual_step, expected_step, tolerance));

  ASSERT_CALL(sleqp_newton_data_release(&newton_data));

  ASSERT_CALL(sleqp_aug_jacobian_release(&jacobian));

  ASSERT_CALL(sleqp_sparse_factorization_release(&factorization));

  ASSERT_CALL(sleqp_options_release(&options));

  ASSERT_CALL(sleqp_params_release(&params));

  ASSERT_CALL(sleqp_sparse_vector_free(&expected_step));

  ASSERT_CALL(sleqp_sparse_vector_free(&actual_step));
}
END_TEST

void newton_teardown()
{
  ASSERT_CALL(sleqp_iterate_release(&iterate));

  ASSERT_CALL(sleqp_problem_free(&problem));

  ASSERT_CALL(sleqp_params_release(&params));

  ASSERT_CALL(sleqp_sparse_vector_free(&linquadfunc_x));

  ASSERT_CALL(sleqp_sparse_vector_free(&linquadfunc_cons_ub));

  ASSERT_CALL(sleqp_sparse_vector_free(&linquadfunc_cons_lb));

  ASSERT_CALL(sleqp_sparse_vector_free(&linquadfunc_var_ub));

  ASSERT_CALL(sleqp_sparse_vector_free(&linquadfunc_var_lb));

  ASSERT_CALL(sleqp_func_release(&linquadfunc));

  sleqp_free(&func_data->x);

  sleqp_free(&func_data);
}

Suite* newton_test_suite()
{
  Suite *suite;
  TCase *tc_cons;

  suite = suite_create("Unconstrained newton step tests");

  tc_cons = tcase_create("Newton step");

  tcase_add_checked_fixture(tc_cons,
                            newton_setup,
                            newton_teardown);

  tcase_add_test(tc_cons, newton_constrained_step);

  suite_add_tcase(suite, tc_cons);

  return suite;
}

int main()
{
  int num_fails;
  Suite* suite;
  SRunner* srunner;

  suite = newton_test_suite();
  srunner = srunner_create(suite);

  srunner_set_fork_status(srunner, CK_NOFORK);
  srunner_run_all(srunner, CK_NORMAL);

  num_fails = srunner_ntests_failed(srunner);

  srunner_free(srunner);

  return (num_fails > 0) ? EXIT_FAILURE : EXIT_SUCCESS;
}
