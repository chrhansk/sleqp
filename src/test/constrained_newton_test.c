#include <check.h>
#include <stdlib.h>

#include "direction.h"
#include "test_common.h"

#include "cmp.h"
#include "func.h"
#include "iterate.h"
#include "mem.h"
#include "newton.h"
#include "problem.h"
#include "util.h"
#include "working_set.h"
#include "working_step.h"

#include "aug_jac/standard_aug_jac.h"
#include "fact/fact.h"

static const int num_variables   = 2;
static const int num_constraints = 1;

SleqpFunc* linquadfunc;

SleqpVec* linquadfunc_var_lb;
SleqpVec* linquadfunc_var_ub;
SleqpVec* linquadfunc_cons_lb;
SleqpVec* linquadfunc_cons_ub;
SleqpVec* linquadfunc_x;

SleqpSettings* settings;
SleqpProblem* problem;
SleqpIterate* iterate;

typedef struct LinQuadFuncData
{
  double* x;
} LinQuadFuncData;

LinQuadFuncData* func_data;

static inline double
square(double v)
{
  return v * v;
}

SLEQP_RETCODE
linquadfunc_set(SleqpFunc* func,
                SleqpVec* x,
                SLEQP_VALUE_REASON reason,
                bool* reject,
                void* func_data)
{
  LinQuadFuncData* data = (LinQuadFuncData*)func_data;

  data->x[0] = 0;
  data->x[1] = 0;

  int k_x = 0;

  while (k_x < x->nnz)
  {
    data->x[x->indices[k_x]] = x->data[k_x];

    ++k_x;
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
linquadfunc_obj_val(SleqpFunc* func, double* obj_val, void* func_data)
{
  LinQuadFuncData* data = (LinQuadFuncData*)func_data;

  *obj_val = square(data->x[0]) + square(data->x[1]);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
linquadfunc_obj_grad(SleqpFunc* func, SleqpVec* obj_grad, void* func_data)
{
  LinQuadFuncData* data = (LinQuadFuncData*)func_data;

  assert(obj_grad->dim == 2);
  assert(obj_grad->nnz_max >= 2);

  obj_grad->nnz = 0;

  SLEQP_CALL(sleqp_vec_push(obj_grad, 0, 2. * data->x[0]));

  SLEQP_CALL(sleqp_vec_push(obj_grad, 1, 2. * data->x[1]));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
linquadfunc_cons_val(SleqpFunc* func, SleqpVec* cons_val, void* func_data)
{
  LinQuadFuncData* data = (LinQuadFuncData*)func_data;

  SLEQP_CALL(sleqp_vec_push(cons_val, 0, data->x[1]));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
linquadfunc_cons_jac(SleqpFunc* func, SleqpMat* cons_jac, void* func_data)
{
  assert(sleqp_mat_nnz(cons_jac) == 0);
  assert(sleqp_mat_nnz_max(cons_jac) >= 1);

  SLEQP_CALL(sleqp_mat_push(cons_jac, 0, 1, 1.));

  assert(sleqp_mat_is_valid(cons_jac));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
linquadfunc_hess_prod(SleqpFunc* func,
                      const SleqpVec* direction,
                      const SleqpVec* cons_duals,
                      SleqpVec* result,
                      void* func_data)
{
  SLEQP_CALL(sleqp_vec_copy(direction, result));

  SLEQP_CALL(sleqp_vec_scale(result, 2.));

  return SLEQP_OKAY;
}

void
newton_setup()
{
  // setup function

  ASSERT_CALL(sleqp_malloc(&func_data));

  ASSERT_CALL(sleqp_alloc_array(&func_data->x, 2));

  SleqpFuncCallbacks callbacks = {.set_value = linquadfunc_set,
                                  .obj_val   = linquadfunc_obj_val,
                                  .obj_grad  = linquadfunc_obj_grad,
                                  .cons_val  = linquadfunc_cons_val,
                                  .cons_jac  = linquadfunc_cons_jac,
                                  .hess_prod = linquadfunc_hess_prod,
                                  .func_free = NULL};

  ASSERT_CALL(sleqp_func_create(&linquadfunc,
                                &callbacks,
                                num_variables,
                                num_constraints,
                                func_data));

  double inf = sleqp_infinity();

  ASSERT_CALL(sleqp_vec_create_full(&linquadfunc_var_lb, num_variables));
  ASSERT_CALL(sleqp_vec_fill(linquadfunc_var_lb, -inf));

  ASSERT_CALL(sleqp_vec_create_full(&linquadfunc_var_ub, num_variables));
  ASSERT_CALL(sleqp_vec_fill(linquadfunc_var_ub, inf));

  ASSERT_CALL(sleqp_vec_create_full(&linquadfunc_cons_lb, num_constraints));

  ASSERT_CALL(sleqp_vec_push(linquadfunc_cons_lb, 0, 2.));

  ASSERT_CALL(sleqp_vec_create_full(&linquadfunc_cons_ub, num_constraints));

  ASSERT_CALL(sleqp_vec_push(linquadfunc_cons_ub, 0, inf));

  ASSERT_CALL(sleqp_vec_create_full(&linquadfunc_x, num_variables));

  ASSERT_CALL(sleqp_vec_push(linquadfunc_x, 0, 1.));

  ASSERT_CALL(sleqp_vec_push(linquadfunc_x, 1, 2.));

  ASSERT_CALL(sleqp_settings_create(&settings));

  ASSERT_CALL(sleqp_problem_create_simple(&problem,
                                          linquadfunc,
                                          settings,
                                          linquadfunc_var_lb,
                                          linquadfunc_var_ub,
                                          linquadfunc_cons_lb,
                                          linquadfunc_cons_ub));

  ASSERT_CALL(sleqp_iterate_create(&iterate, problem, linquadfunc_x));

  ASSERT_CALL(
    sleqp_set_and_evaluate(problem, iterate, SLEQP_VALUE_REASON_NONE, NULL));

  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);
  ASSERT_CALL(sleqp_working_set_reset(working_set));

  // set the cons state manually...

  ASSERT_CALL(sleqp_working_set_add_cons(working_set, 0, SLEQP_ACTIVE_LOWER));
}

START_TEST(newton_constrained_step)
{
  SleqpVec* expected_step;

  SleqpDirection* actual_direction;

  SleqpSettings* settings;

  SleqpWorkingStep* working_step;
  SleqpEQPSolver* newton_solver;

  SleqpFact* fact;
  SleqpAugJac* jacobian;

  double penalty_parameter = 1.;
  double trust_radius      = 10.;

  const int num_variables = sleqp_problem_num_vars(problem);

  ASSERT_CALL(sleqp_vec_create(&expected_step, num_variables, 1));

  ASSERT_CALL(sleqp_vec_push(expected_step, 0, -1.));

  ASSERT_CALL(sleqp_settings_create(&settings));

  ASSERT_CALL(sleqp_direction_create(&actual_direction, problem, settings));

  SleqpVec* actual_step = sleqp_direction_primal(actual_direction);

  ASSERT_CALL(sleqp_fact_create_default(&fact, settings));

  ASSERT_CALL(sleqp_standard_aug_jac_create(&jacobian, problem, settings, fact));

  ASSERT_CALL(sleqp_aug_jac_set_iterate(jacobian, iterate));

  ASSERT_CALL(sleqp_working_step_create(&working_step, problem, settings));

  ASSERT_CALL(sleqp_newton_solver_create(&newton_solver,
                                         problem,
                                         settings,
                                         working_step));

  ASSERT_CALL(sleqp_eqp_solver_set_iterate(newton_solver,
                                           iterate,
                                           jacobian,
                                           trust_radius,
                                           penalty_parameter));

  ASSERT_CALL(
    sleqp_eqp_solver_compute_direction(newton_solver,
                                       sleqp_iterate_cons_dual(iterate),
                                       actual_direction));

  const double tolerance = 1e-8;

  ck_assert(sleqp_vec_eq(actual_step, expected_step, tolerance));

  ASSERT_CALL(sleqp_eqp_solver_release(&newton_solver));

  ASSERT_CALL(sleqp_working_step_release(&working_step));

  ASSERT_CALL(sleqp_aug_jac_release(&jacobian));

  ASSERT_CALL(sleqp_fact_release(&fact));

  ASSERT_CALL(sleqp_direction_release(&actual_direction));

  ASSERT_CALL(sleqp_settings_release(&settings));

  ASSERT_CALL(sleqp_vec_free(&expected_step));
}
END_TEST

void
newton_teardown()
{
  ASSERT_CALL(sleqp_iterate_release(&iterate));

  ASSERT_CALL(sleqp_problem_release(&problem));

  ASSERT_CALL(sleqp_settings_release(&settings));

  ASSERT_CALL(sleqp_vec_free(&linquadfunc_x));

  ASSERT_CALL(sleqp_vec_free(&linquadfunc_cons_ub));

  ASSERT_CALL(sleqp_vec_free(&linquadfunc_cons_lb));

  ASSERT_CALL(sleqp_vec_free(&linquadfunc_var_ub));

  ASSERT_CALL(sleqp_vec_free(&linquadfunc_var_lb));

  ASSERT_CALL(sleqp_func_release(&linquadfunc));

  sleqp_free(&func_data->x);

  sleqp_free(&func_data);
}

Suite*
newton_test_suite()
{
  Suite* suite;
  TCase* tc_cons;

  suite = suite_create("Unconstrained newton step tests");

  tc_cons = tcase_create("Newton step");

  tcase_add_checked_fixture(tc_cons, newton_setup, newton_teardown);

  tcase_add_test(tc_cons, newton_constrained_step);

  suite_add_tcase(suite, tc_cons);

  return suite;
}

TEST_MAIN(newton_test_suite)
