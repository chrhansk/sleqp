#include <stdlib.h>
#include <check.h>

#include "sleqp.h"
#include "sleqp_cauchy.h"
#include "sleqp_cmp.h"
#include "sleqp_mem.h"

#include "lp/sleqp_lpi_soplex.h"

#include "test_common.h"

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

static SLEQP_RETCODE linfunc_set(SleqpSparseVec* x,
                                 int num_variables,
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

static SLEQP_RETCODE linfunc_eval(int num_variables,
                                  SleqpSparseVec* indices,
                                  double* func_val,
                                  SleqpSparseVec* func_grad,
                                  SleqpSparseVec* cons_val,
                                  SleqpSparseMatrix* cons_jac,
                                  void* func_data)
{
  LinFuncData* data = (LinFuncData*) func_data;

  if(func_val)
  {
    *func_val = data->x[0] + data->x[1];
  }

  if(func_grad)
  {
    assert(func_grad->dim == 2);

    func_grad->nnz = 0;

    SLEQP_CALL(sleqp_sparse_vector_push(func_grad,
                                        0,
                                        1.));

    SLEQP_CALL(sleqp_sparse_vector_push(func_grad,
                                        1,
                                        1.));

  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE linfunc_hess_prod(int num_variables,
                                       double* func_dual,
                                       SleqpSparseVec* direction,
                                       SleqpSparseVec* cons_duals,
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

  ASSERT_CALL(sleqp_calloc(&func_data->x, 2));

  ASSERT_CALL(sleqp_func_create(&linfunc,
                                linfunc_set,
                                linfunc_eval,
                                linfunc_hess_prod,
                                2,
                                func_data));

  ASSERT_CALL(sleqp_sparse_vector_create(&linfunc_var_lb,
                                         2,
                                         2));

  ASSERT_CALL(sleqp_sparse_vector_push(linfunc_var_lb, 0, -inf));
  ASSERT_CALL(sleqp_sparse_vector_push(linfunc_var_lb, 1, -inf));

  ASSERT_CALL(sleqp_sparse_vector_create(&linfunc_var_ub,
                                         2,
                                         2));

  ASSERT_CALL(sleqp_sparse_vector_push(linfunc_var_ub, 0, inf));
  ASSERT_CALL(sleqp_sparse_vector_push(linfunc_var_ub, 1, inf));

  ASSERT_CALL(sleqp_sparse_vector_create(&linfunc_cons_lb,
                                         0,
                                         0));

  ASSERT_CALL(sleqp_sparse_vector_create(&linfunc_cons_ub,
                                         0,
                                         0));

  ASSERT_CALL(sleqp_sparse_vector_create(&linfunc_x,
                                         2,
                                         2));

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



  ASSERT_CALL(sleqp_func_free(&linfunc));

  sleqp_free(&func_data->x);

  sleqp_free(&func_data);
}


START_TEST(test_unconstrained_cauchy_direction)
{
  SleqpProblem* problem;
  SleqpIterate* iterate;
  SleqpLPi* lp_interface;
  SleqpSparseVec* direction;
  SleqpCauchyData* cauchy_data;

  double penalty = 1., trust_radius = 1.5;


  ASSERT_CALL(sleqp_problem_create(&problem,
                                   linfunc,
                                   linfunc_var_lb,
                                   linfunc_var_ub,
                                   linfunc_cons_lb,
                                   linfunc_cons_ub));

  ASSERT_CALL(sleqp_iterate_create(&iterate,
                                   problem,
                                   linfunc_x));

  int num_lp_variables = problem->num_variables + 2*problem->num_constraints;
  int num_lp_constraints = problem->num_constraints;

  ASSERT_CALL(sleqp_lpi_soplex_create_interface(&lp_interface,
                                                num_lp_variables,
                                                num_lp_constraints));

  ASSERT_CALL(sleqp_set_and_evaluate(problem, iterate));

  ASSERT_CALL(sleqp_sparse_vector_create(&direction, 0, 0));

  ASSERT_CALL(sleqp_cauchy_data_create(&cauchy_data,
                                       problem,
                                       lp_interface));

  ASSERT_CALL(sleqp_cauchy_compute_direction(cauchy_data,
                                             iterate,
                                             penalty,
                                             trust_radius));

  ASSERT_CALL(sleqp_cauchy_get_direction(cauchy_data,
                                         iterate,
                                         direction));

  ck_assert_int_eq(direction->dim, 2);

  ck_assert(sleqp_sparse_vector_at(direction, 0));
  ck_assert(sleqp_sparse_vector_at(direction, 1));

  ck_assert(sleqp_eq(*sleqp_sparse_vector_at(direction, 0), -trust_radius));
  ck_assert(sleqp_eq(*sleqp_sparse_vector_at(direction, 1), -trust_radius));

  ASSERT_CALL(sleqp_cauchy_data_free(&cauchy_data));

  ASSERT_CALL(sleqp_sparse_vector_free(&direction));

  ASSERT_CALL(sleqp_lpi_free(&lp_interface));

  ASSERT_CALL(sleqp_iterate_free(&iterate));

  ASSERT_CALL(sleqp_problem_free(&problem));
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


int main()
{
  int num_fails;
  Suite* suite;
  SRunner* srunner;

  suite = cauchy_test_suite();
  srunner = srunner_create(suite);

  srunner_set_fork_status(srunner, CK_NOFORK);
  srunner_run_all(srunner, CK_NORMAL);

  num_fails = srunner_ntests_failed(srunner);

  srunner_free(srunner);

  return (num_fails > 0) ? EXIT_FAILURE : EXIT_SUCCESS;
}
