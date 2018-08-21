#include <stdlib.h>
#include <check.h>

#include "sleqp.h"
#include "sleqp_cauchy.h"
#include "sleqp_cmp.h"
#include "sleqp_dual_estimation.h"
#include "sleqp_mem.h"

#include "lp/sleqp_lpi_soplex.h"

#include "test_common.h"

SleqpFunc* quadfunc;

SleqpSparseVec* var_lb;
SleqpSparseVec* var_ub;
SleqpSparseVec* cons_lb;
SleqpSparseVec* cons_ub;
SleqpSparseVec* x;

typedef struct SquareFuncData
{
  double* x;
} SquareFuncData;

SquareFuncData* func_data;

SleqpSparseVec* quadfunc_var_lb;
SleqpSparseVec* quadfunc_var_ub;
SleqpSparseVec* quadfunc_cons_lb;
SleqpSparseVec* quadfunc_cons_ub;
SleqpSparseVec* quadfunc_x;

static inline double square(double v)
{
  return v*v;
}

static SLEQP_RETCODE quadfunc_set(SleqpSparseVec* x,
                                 size_t num_variables,
                                 size_t* func_grad_nnz,
                                 size_t* cons_val_nnz,
                                 size_t* cons_jac_nnz,
                                 void* func_data)
{
  *func_grad_nnz = 2;
  *cons_val_nnz = 0;
  *cons_jac_nnz = 0;

  SquareFuncData* data = (SquareFuncData*) func_data;

  data->x[0] = 0;
  data->x[1] = 0;

  size_t k_x = 0;

  while(k_x < x->nnz)
  {
    data->x[x->indices[k_x]] = x->data[k_x];

    ++k_x;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE quadfunc_eval(size_t num_variables,
                                  int* indices,
                                  double* func_val,
                                  SleqpSparseVec* func_grad,
                                  SleqpSparseVec* cons_val,
                                  SleqpSparseMatrix* cons_jac,
                                  void* func_data)
{
  SquareFuncData* data = (SquareFuncData*) func_data;

  if(func_val)
  {
    *func_val = square(data->x[0]) + square(data->x[1]);
  }

  if(func_grad)
  {
    assert(func_grad->dim == 2);

    func_grad->nnz = 0;

    SLEQP_CALL(sleqp_sparse_vector_push(func_grad,
                                        0,
                                        2.*data->x[0]));

    SLEQP_CALL(sleqp_sparse_vector_push(func_grad,
                                        1,
                                        2.*data->x[1]));

  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE quadfunc_eval_bilinear(size_t num_variables,
                                           double* fval,
                                           SleqpSparseVec* direction,
                                           SleqpSparseVec* multipliers,
                                           void* func_data)
{
  return SLEQP_OKAY;
}

void simply_constrained_setup()
{
  const double inf = sleqp_infinity();

  ASSERT_CALL(sleqp_malloc(&func_data));

  ASSERT_CALL(sleqp_calloc(&func_data->x, 2));

  ASSERT_CALL(sleqp_func_create(&quadfunc,
                                quadfunc_set,
                                quadfunc_eval,
                                quadfunc_eval_bilinear,
                                func_data));

  ASSERT_CALL(sleqp_sparse_vector_create(&quadfunc_var_lb,
                                         2,
                                         2));

  ASSERT_CALL(sleqp_sparse_vector_push(quadfunc_var_lb, 0, 1));
  ASSERT_CALL(sleqp_sparse_vector_push(quadfunc_var_lb, 1, 2));

  ASSERT_CALL(sleqp_sparse_vector_create(&quadfunc_var_ub,
                                         2,
                                         2));

  ASSERT_CALL(sleqp_sparse_vector_push(quadfunc_var_ub, 0, 2));
  ASSERT_CALL(sleqp_sparse_vector_push(quadfunc_var_ub, 1, 3));

  ASSERT_CALL(sleqp_sparse_vector_create(&quadfunc_cons_lb,
                                         0,
                                         0));

  ASSERT_CALL(sleqp_sparse_vector_create(&quadfunc_cons_ub,
                                         0,
                                         0));

  ASSERT_CALL(sleqp_sparse_vector_create(&quadfunc_x,
                                         2,
                                         2));

  ASSERT_CALL(sleqp_sparse_vector_push(quadfunc_x, 0, 1));
  ASSERT_CALL(sleqp_sparse_vector_push(quadfunc_x, 1, 2));
}

void simply_constrained_teardown()
{
  ASSERT_CALL(sleqp_sparse_vector_free(&quadfunc_x));

  ASSERT_CALL(sleqp_sparse_vector_free(&quadfunc_cons_ub));

  ASSERT_CALL(sleqp_sparse_vector_free(&quadfunc_cons_lb));

  ASSERT_CALL(sleqp_sparse_vector_free(&quadfunc_var_ub));

  ASSERT_CALL(sleqp_sparse_vector_free(&quadfunc_var_lb));


  ASSERT_CALL(sleqp_func_free(&quadfunc));

  sleqp_free(&func_data->x);

  sleqp_free(&func_data);
}

START_TEST(test_simply_constrained_dual_estimation)
{
  SleqpProblem* problem;
  SleqpIterate* iterate;
  SleqpLPi* lp_interface;
  SleqpCauchyData* cauchy_data;
  SleqpActiveSet* active_set;
  SleqpDualEstimationData* estimation_data;

  size_t num_variables;
  size_t num_constraints;

  double penalty = 1., trust_radius = 1.5;


  ASSERT_CALL(sleqp_problem_create(&problem,
                                   quadfunc,
                                   quadfunc_var_lb,
                                   quadfunc_var_ub,
                                   quadfunc_cons_lb,
                                   quadfunc_cons_ub));

  num_variables = problem->num_variables;
  num_constraints = problem->num_constraints;

  ASSERT_CALL(sleqp_iterate_create(&iterate,
                                   problem,
                                   quadfunc_x));

  ASSERT_CALL(sleqp_lpi_soplex_create_interface(&lp_interface));

  ASSERT_CALL(sleqp_lpi_create_problem(lp_interface,
                                       num_variables + 2*num_constraints,
                                       num_constraints));

  ASSERT_CALL(sleqp_set_and_evaluate(problem, iterate));

  ASSERT_CALL(sleqp_cauchy_data_create(&cauchy_data,
                                       problem,
                                       lp_interface));

  ASSERT_CALL(sleqp_active_set_create(&active_set,
                                      problem));

  ASSERT_CALL(sleqp_dual_estimation_data_create(&estimation_data, problem));

  ASSERT_CALL(sleqp_cauchy_compute_direction(problem,
                                             iterate,
                                             cauchy_data,
                                             penalty,
                                             trust_radius));

  ASSERT_CALL(sleqp_cauchy_get_active_set(problem,
                                          iterate,
                                          cauchy_data,
                                          active_set,
                                          trust_radius));

  ASSERT_CALL(sleqp_dual_estimation_compute(estimation_data,
                                            iterate,
                                            active_set));

  SleqpSparseVec* vars_dual = iterate->vars_dual;

  ck_assert(sleqp_sparse_vector_at(vars_dual, 0));
  ck_assert(sleqp_sparse_vector_at(vars_dual, 1));

  ck_assert(sleqp_eq(*sleqp_sparse_vector_at(vars_dual, 0), -2.));
  ck_assert(sleqp_eq(*sleqp_sparse_vector_at(vars_dual, 1), -4.));

  ASSERT_CALL(sleqp_dual_estimation_data_free(&estimation_data));

  ASSERT_CALL(sleqp_active_set_free(&active_set));

  ASSERT_CALL(sleqp_cauchy_data_free(&cauchy_data));

  ASSERT_CALL(sleqp_lpi_free(&lp_interface));

  ASSERT_CALL(sleqp_iterate_free(&iterate));

  ASSERT_CALL(sleqp_problem_free(&problem));
}
END_TEST

Suite* dual_estimation_test_suite()
{
  Suite *suite;
  TCase *tc_dual_estimation;

  suite = suite_create("Dual estimation tests");

  tc_dual_estimation = tcase_create("Simply constrained");

  tcase_add_checked_fixture(tc_dual_estimation,
                            simply_constrained_setup,
                            simply_constrained_teardown);

  tcase_add_test(tc_dual_estimation, test_simply_constrained_dual_estimation);
  suite_add_tcase(suite, tc_dual_estimation);

  return suite;
}

int main()
{
  int num_fails;
  Suite* suite;
  SRunner* srunner;

  suite = dual_estimation_test_suite();
  srunner = srunner_create(suite);

  srunner_set_fork_status(srunner, CK_NOFORK);
  srunner_run_all(srunner, CK_NORMAL);

  num_fails = srunner_ntests_failed(srunner);

  srunner_free(srunner);

  return (num_fails > 0) ? EXIT_FAILURE : EXIT_SUCCESS;
}
