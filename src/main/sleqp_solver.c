#include "sleqp_solver.h"

#include <assert.h>

#include "sleqp_cmp.h"
#include "sleqp_mem.h"

#include "sleqp_defs.h"

#include "sleqp_aug_jacobian.h"

#include "sleqp_cauchy.h"
#include "sleqp_dual_estimation.h"

#include "sleqp_newton.h"

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

  SleqpSparseVec* cauchy_hessian_direction;

  SleqpSparseVec* cauchy_step;

  SleqpNewtonData* newton_data;

  SleqpSparseVec* newton_step;

  SleqpSparseVec* cauchy_newton_direction;

  SleqpSparseVec* cauchy_newton_hessian_direction;

  SleqpSparseVec* trial_point;

  SleqpAugJacobian* aug_jacobian;

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

  assert(sleqp_sparse_vector_valid(x));

  SLEQP_CALL(sleqp_sparse_vector_clip(x,
                                      solver->problem->var_lb,
                                      solver->problem->var_ub,
                                      &xclip));

  SLEQP_CALL(sleqp_iterate_create(&solver->iterate,
                                  solver->problem,
                                  xclip));

  // TODO: make this generic at a later point...

  int num_lp_variables = num_variables + 2*num_constraints;
  int num_lp_constraints = num_constraints;

  SLEQP_CALL(sleqp_lpi_soplex_create_interface(&solver->lp_interface,
                                               num_lp_variables,
                                               num_lp_constraints));

  SLEQP_CALL(sleqp_cauchy_data_create(&solver->cauchy_data,
                                      problem,
                                      solver->lp_interface));

  SLEQP_CALL(sleqp_sparse_vector_create(&solver->cauchy_direction,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&solver->cauchy_hessian_direction,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&solver->cauchy_step,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_newton_data_create(&solver->newton_data,
                                      problem));

  SLEQP_CALL(sleqp_sparse_vector_create(&solver->newton_step,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&solver->cauchy_newton_direction,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&solver->cauchy_newton_hessian_direction,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&solver->trial_point,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_aug_jacobian_create(&solver->aug_jacobian, problem));

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

static SLEQP_RETCODE compute_trial_point(SleqpSolver* solver)
{
  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;

  double one = 1.;

  SLEQP_CALL(sleqp_set_and_evaluate(problem, iterate));

  // compute Cauchy direction / step and dual estimation
  {
  SLEQP_CALL(sleqp_cauchy_compute_direction(solver->cauchy_data,
                                            iterate,
                                            solver->penalty_parameter,
                                            solver->trust_radius));

  SLEQP_CALL(sleqp_cauchy_get_active_set(solver->cauchy_data,
                                         iterate,
                                         solver->trust_radius));

  SLEQP_CALL(sleqp_aug_jacobian_set_iterate(solver->aug_jacobian,
                                            iterate));

  SLEQP_CALL(sleqp_cauchy_get_direction(solver->cauchy_data,
                                        iterate,
                                        solver->cauchy_direction));

  SLEQP_CALL(sleqp_dual_estimation_compute(solver->estimation_data,
                                           iterate,
                                           solver->aug_jacobian));

  SLEQP_CALL(sleqp_func_hess_product(problem->func,
                                     &one,
                                     solver->cauchy_direction,
                                     iterate->cons_dual,
                                     solver->cauchy_hessian_direction));

  SLEQP_CALL(sleqp_cauchy_compute_step(solver->cauchy_data,
                                       iterate,
                                       solver->penalty_parameter,
                                       solver->cauchy_hessian_direction,
                                       solver->cauchy_step));
  }

  // compute Newton step
  {
    SLEQP_CALL(sleqp_newton_compute_step(solver->newton_data,
                                         solver->iterate,
                                         solver->aug_jacobian,
                                         solver->trust_radius,
                                         solver->penalty_parameter,
                                         solver->newton_step));

  }

  // compute Cauchy Newton direction (from Cauchy point to Newton point)
  {
    SLEQP_CALL(sleqp_sparse_vector_add(solver->newton_step,
                                       solver->cauchy_step,
                                       1.,
                                       -1.,
                                       solver->cauchy_newton_direction));

    SLEQP_CALL(sleqp_func_hess_product(problem->func,
                                       &one,
                                       solver->cauchy_newton_direction,
                                       iterate->cons_dual,
                                       solver->cauchy_newton_hessian_direction));
  }


  return SLEQP_OKAY;
}


SLEQP_RETCODE sleqp_solve(SleqpSolver* solver)
{
  SLEQP_CALL(compute_trial_point(solver));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solver_free(SleqpSolver** star)
{
  SleqpSolver* solver = *star;

  SLEQP_CALL(sleqp_dual_estimation_data_free(&solver->estimation_data));

  SLEQP_CALL(sleqp_aug_jacobian_free(&solver->aug_jacobian));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->trial_point));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->cauchy_newton_hessian_direction));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->cauchy_newton_direction));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->newton_step));

  SLEQP_CALL(sleqp_newton_data_free(&solver->newton_data));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->cauchy_step));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->cauchy_hessian_direction));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->cauchy_direction));

  SLEQP_CALL(sleqp_cauchy_data_free(&solver->cauchy_data));

  SLEQP_CALL(sleqp_lpi_free(&solver->lp_interface));

  SLEQP_CALL(sleqp_iterate_free(&solver->iterate));

  sleqp_free(star);

  return SLEQP_OKAY;
}
