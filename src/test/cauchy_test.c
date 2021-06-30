#include <stdlib.h>
#include <check.h>

#include "sleqp.h"
#include "sleqp_cauchy.h"
#include "sleqp_cmp.h"
#include "sleqp_mem.h"

#include "lp/sleqp_lpi.h"

#include "test_common.h"

const int num_variables = 2;
const int num_constraints = 0;

SleqpFunc* linfunc;

SleqpSparseVec* var_lb;
SleqpSparseVec* var_ub;
SleqpSparseVec* cons_lb;
SleqpSparseVec* cons_ub;
SleqpSparseVec* x;

typedef struct LinFuncData
{
  double* x;
} LinFuncData;

LinFuncData* func_data;

SleqpSparseVec* linfunc_var_lb;
SleqpSparseVec* linfunc_var_ub;
SleqpSparseVec* linfunc_cons_lb;
SleqpSparseVec* linfunc_cons_ub;
SleqpSparseVec* linfunc_x;

static SLEQP_RETCODE linfunc_set(SleqpFunc* func,
                                 SleqpSparseVec* x,
                                 SLEQP_VALUE_REASON reason,
                                 int* func_grad_nnz,
                                 int* cons_val_nnz,
                                 int* cons_jac_nnz,
                                 void* func_data)
{
  *func_grad_nnz = 2;
  *cons_val_nnz = 0;
  *cons_jac_nnz = 0;

  LinFuncData* data = (LinFuncData*) func_data;

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

static SLEQP_RETCODE linfunc_val(SleqpFunc* func,
                                 double* func_val,
                                 void* func_data)
{
  LinFuncData* data = (LinFuncData*) func_data;

  *func_val = data->x[0] + data->x[1];

  return SLEQP_OKAY;
}

static SLEQP_RETCODE linfunc_grad(SleqpFunc* func,
                                  SleqpSparseVec* func_grad,
                                  void* func_data)
{
  assert(func_grad->dim == 2);

  func_grad->nnz = 0;

  SLEQP_CALL(sleqp_sparse_vector_push(func_grad,
                                      0,
                                      1.));

  SLEQP_CALL(sleqp_sparse_vector_push(func_grad,
                                      1,
                                      1.));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE linfunc_hess_prod(SleqpFunc* func,
                                       const double* func_dual,
                                       const SleqpSparseVec* direction,
                                       const SleqpSparseVec* cons_duals,
                                       SleqpSparseVec* product,
                                       void* func_data)
{
  product->nnz = 0;
  return SLEQP_OKAY;
}


void unconstrained_setup()
{
  const double inf = sleqp_infinity();

  ASSERT_CALL(sleqp_malloc(&func_data));

  ASSERT_CALL(sleqp_alloc_array(&func_data->x, 2));

  SleqpFuncCallbacks callbacks = {
    .set_value = linfunc_set,
    .func_val  = linfunc_val,
    .func_grad = linfunc_grad,
    .cons_val  = NULL,
    .cons_jac  = NULL,
    .hess_prod = linfunc_hess_prod,
    .func_free = NULL
  };

  ASSERT_CALL(sleqp_func_create(&linfunc,
                                &callbacks,
                                num_variables,
                                num_constraints,
                                func_data));

  ASSERT_CALL(sleqp_sparse_vector_create_full(&linfunc_var_lb,
                                              num_variables));

  ASSERT_CALL(sleqp_sparse_vector_push(linfunc_var_lb, 0, -inf));
  ASSERT_CALL(sleqp_sparse_vector_push(linfunc_var_lb, 1, -inf));

  ASSERT_CALL(sleqp_sparse_vector_create_full(&linfunc_var_ub,
                                              num_variables));

  ASSERT_CALL(sleqp_sparse_vector_push(linfunc_var_ub, 0, inf));
  ASSERT_CALL(sleqp_sparse_vector_push(linfunc_var_ub, 1, inf));

  ASSERT_CALL(sleqp_sparse_vector_create_empty(&linfunc_cons_lb,
                                               num_constraints));

  ASSERT_CALL(sleqp_sparse_vector_create_empty(&linfunc_cons_ub,
                                               num_constraints));

  ASSERT_CALL(sleqp_sparse_vector_create_full(&linfunc_x,
                                              num_variables));

  ASSERT_CALL(sleqp_sparse_vector_push(linfunc_x, 0, 0.));
  ASSERT_CALL(sleqp_sparse_vector_push(linfunc_x, 1, 0.));
}

void unconstrained_teardown()
{
  ASSERT_CALL(sleqp_sparse_vector_free(&linfunc_x));

  ASSERT_CALL(sleqp_sparse_vector_free(&linfunc_cons_ub));

  ASSERT_CALL(sleqp_sparse_vector_free(&linfunc_cons_lb));

  ASSERT_CALL(sleqp_sparse_vector_free(&linfunc_var_ub));

  ASSERT_CALL(sleqp_sparse_vector_free(&linfunc_var_lb));



  ASSERT_CALL(sleqp_func_release(&linfunc));

  sleqp_free(&func_data->x);

  sleqp_free(&func_data);
}


START_TEST(test_unconstrained_cauchy_direction)
{
  SleqpParams* params;
  SleqpOptions* options;
  SleqpProblem* problem;
  SleqpIterate* iterate;
  SleqpLPi* lp_interface;
  SleqpSparseVec* direction;
  SleqpCauchy* cauchy_data;

  double penalty_parameter = 1., trust_radius = 1.5;

  ASSERT_CALL(sleqp_params_create(&params));

  ASSERT_CALL(sleqp_options_create(&options));

  ASSERT_CALL(sleqp_problem_create_simple(&problem,
                                          linfunc,
                                          params,
                                          linfunc_var_lb,
                                          linfunc_var_ub,
                                          linfunc_cons_lb,
                                          linfunc_cons_ub));

  const int num_variables = sleqp_problem_num_variables(problem);
  const int num_constraints = sleqp_problem_num_constraints(problem);

  ASSERT_CALL(sleqp_iterate_create(&iterate,
                                   problem,
                                   linfunc_x));

  int num_lp_variables = num_variables + 2*num_constraints;
  int num_lp_constraints = num_constraints;

  ASSERT_CALL(sleqp_lpi_create_default_interface(&lp_interface,
                                                 num_lp_variables,
                                                 num_lp_constraints,
                                                 params));

  ASSERT_CALL(sleqp_set_and_evaluate(problem, iterate, SLEQP_VALUE_REASON_NONE));

  ASSERT_CALL(sleqp_sparse_vector_create(&direction, 0, 0));

  ASSERT_CALL(sleqp_cauchy_create(&cauchy_data,
                                  problem,
                                  params,
                                  options,
                                  lp_interface));

  ASSERT_CALL(sleqp_cauchy_set_iterate(cauchy_data,
                                       iterate,
                                       trust_radius));

  ASSERT_CALL(sleqp_cauchy_solve(cauchy_data,
                                 sleqp_iterate_get_func_grad(iterate),
                                 penalty_parameter,
                                 SLEQP_CAUCHY_OBJECTIVE_TYPE_DEFAULT));

  ASSERT_CALL(sleqp_cauchy_get_direction(cauchy_data,
                                         direction));

  ck_assert_int_eq(direction->dim, 2);

  ck_assert(sleqp_sparse_vector_at(direction, 0));
  ck_assert(sleqp_sparse_vector_at(direction, 1));

  double tolerance = 1e-8;

  ck_assert(sleqp_is_eq(*sleqp_sparse_vector_at(direction, 0), -trust_radius, tolerance));
  ck_assert(sleqp_is_eq(*sleqp_sparse_vector_at(direction, 1), -trust_radius, tolerance));

  ASSERT_CALL(sleqp_cauchy_free(&cauchy_data));

  ASSERT_CALL(sleqp_sparse_vector_free(&direction));

  ASSERT_CALL(sleqp_lpi_free(&lp_interface));

  ASSERT_CALL(sleqp_iterate_release(&iterate));

  ASSERT_CALL(sleqp_problem_release(&problem));

  ASSERT_CALL(sleqp_options_release(&options));

  ASSERT_CALL(sleqp_params_release(&params));
}
END_TEST

Suite* cauchy_test_suite()
{
  Suite *suite;
  TCase *tc_uncons;

  suite = suite_create("Cauchy tests");

  tc_uncons = tcase_create("Unconstrained Cauchy direction");

  tcase_add_checked_fixture(tc_uncons,
                            unconstrained_setup,
                            unconstrained_teardown);

  tcase_add_test(tc_uncons, test_unconstrained_cauchy_direction);
  suite_add_tcase(suite, tc_uncons);

  return suite;
}

TEST_MAIN(cauchy_test_suite)
