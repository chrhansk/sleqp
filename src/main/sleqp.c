#include "sleqp.h"

#include <assert.h>

#include "sleqp_cmp.h"
#include "sleqp_mem.h"

#include "sleqp_cauchy.h"
#include "sleqp_dual_estimation.h"
#include "sleqp_iterate.h"

#include "lp/sleqp_lpi.h"
#include "lp/sleqp_lpi_soplex.h"

struct SleqpSolver
{
  SleqpProblem* problem;

  SleqpIterate* iterate;

  SleqpLPi* lp_interface;

  SleqpCauchyData* cauchy_data;
  SleqpSparseVec* cauchy_direction;

  SleqpSparseVec* trial_point;

  SleqpDualEstimationData* estimation_data;

  double trust_radius;
  double penalty_parameter;
};

SLEQP_RETCODE sleqp_solver_create(SleqpSolver** star,
                                  SleqpProblem* problem,
                                  SleqpSparseVec* x)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpSolver* solver = *star;
  int num_constraints = problem->num_constraints;
  int num_variables = problem->num_variables;

  solver->problem = problem;

  SleqpSparseVec* xclip;

  SLEQP_CALL(sleqp_sparse_vector_clip(x,
                                      solver->problem->var_lb,
                                      solver->problem->var_ub,
                                      &xclip));

  SLEQP_CALL(sleqp_iterate_create(&solver->iterate,
                                  solver->problem,
                                  xclip));

  // TODO: make this generic at a later point...
  SLEQP_CALL(sleqp_lpi_soplex_create_interface(&solver->lp_interface));

  SLEQP_CALL(sleqp_lpi_create_problem(solver->lp_interface,
                                      num_variables + 2*num_constraints,
                                      num_constraints));

  SLEQP_CALL(sleqp_cauchy_data_create(&solver->cauchy_data,
                                      problem,
                                      solver->lp_interface));

  SLEQP_CALL(sleqp_sparse_vector_create(&solver->cauchy_direction,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&solver->trial_point,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_dual_estimation_data_create(&solver->estimation_data,
                                               problem));

  // TODO: Set this to something different?!
  solver->trust_radius = 1.;
  solver->penalty_parameter = 1.;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE update_trust_radius(double reduction_ratio,
                                         double direction_norm,
                                         double* trust_radius)
{
  if(reduction_ratio >= 0.9)
  {
    *trust_radius = SLEQP_MAX(*trust_radius, 7*direction_norm);
  }
  else if(reduction_ratio >= 0.3)
  {
    *trust_radius = SLEQP_MAX(*trust_radius, 2*direction_norm);
  }
  else if(reduction_ratio >= 1e-8)
  {
  }
  else
  {
    *trust_radius = SLEQP_MIN(0.5*(*trust_radius),
                              0.5*direction_norm);
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_set_and_evaluate(SleqpProblem* problem,
                                     SleqpIterate* iterate)
{
  int func_grad_nnz = 0;
  int cons_val_nnz = 0;
  int cons_jac_nnz = 0;

  SLEQP_CALL(sleqp_func_set_value(problem->func,
                                  iterate->x,
                                  &func_grad_nnz,
                                  &cons_val_nnz,
                                  &cons_jac_nnz));

  SLEQP_CALL(sleqp_sparse_vector_reserve(iterate->func_grad, func_grad_nnz));

  SLEQP_CALL(sleqp_sparse_vector_reserve(iterate->cons_val, cons_val_nnz));

  SLEQP_CALL(sleqp_sparse_matrix_reserve(iterate->cons_jac, cons_jac_nnz));

  SLEQP_CALL(sleqp_func_eval(problem->func,
                             NULL,
                             &iterate->func_val,
                             iterate->func_grad,
                             iterate->cons_val,
                             iterate->cons_jac));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE compute_trial_point(SleqpSolver* solver,
                                         double* predicted_reduction)
{
  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;

  SLEQP_CALL(sleqp_set_and_evaluate(problem, iterate));


  SLEQP_CALL(sleqp_cauchy_compute_direction(solver->cauchy_data,
                                            iterate,
                                            solver->penalty_parameter,
                                            solver->trust_radius));

  SLEQP_CALL(sleqp_cauchy_get_active_set(solver->cauchy_data,
                                         iterate,
                                         solver->trust_radius));

  SLEQP_CALL(sleqp_cauchy_get_direction(solver->cauchy_data,
                                        iterate,
                                        solver->cauchy_direction));

  SLEQP_CALL(sleqp_dual_estimation_compute(solver->estimation_data,
                                           iterate));

  SLEQP_CALL(sleqp_cauchy_compute_step(solver->cauchy_data,
                                       iterate,
                                       solver->penalty_parameter,
                                       solver->cauchy_direction,
                                       predicted_reduction));

  SLEQP_CALL(sleqp_sparse_vector_add(iterate->x,
                                     solver->cauchy_direction,
                                     1.,
                                     1.,
                                     solver->trial_point));


  return SLEQP_OKAY;
}


SLEQP_RETCODE sleqp_solve(SleqpSolver* solver)
{
  double predicted_reduction;

  SLEQP_CALL(compute_trial_point(solver, &predicted_reduction));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solver_free(SleqpSolver** star)
{
  SleqpSolver* solver = *star;

  SLEQP_CALL(sleqp_dual_estimation_data_free(&solver->estimation_data));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->trial_point));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->cauchy_direction));

  SLEQP_CALL(sleqp_cauchy_data_free(&solver->cauchy_data));

  SLEQP_CALL(sleqp_lpi_free(&solver->lp_interface));

  SLEQP_CALL(sleqp_iterate_free(&solver->iterate));

  sleqp_free(star);

  return SLEQP_OKAY;
}
