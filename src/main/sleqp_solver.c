#include "sleqp_solver.h"

#include <assert.h>
#include <math.h>

#include "sleqp_cmp.h"
#include "sleqp_mem.h"

#include "sleqp_defs.h"

#include "sleqp_deriv_check.h"

#include "sleqp_aug_jacobian.h"

#include "sleqp_cauchy.h"
#include "sleqp_dual_estimation.h"

#include "sleqp_newton.h"

#include "sleqp_soc.h"

#include "sleqp_iterate.h"
#include "sleqp_merit.h"

#include "lp/sleqp_lpi.h"
#include "lp/sleqp_lpi_soplex.h"

struct SleqpSolver
{
  SleqpProblem* problem;

  SLEQP_STATUS status;

  SleqpParams* params;

  SleqpDerivCheckData* deriv_check;

  SleqpIterate* iterate;

  SleqpSparseVec* violation;

  SleqpLPi* lp_interface;

  SleqpCauchyData* cauchy_data;

  SleqpSparseVec* cauchy_direction;

  SleqpSparseVec* cauchy_hessian_direction;

  SleqpSparseVec* cauchy_step;
  double cauchy_step_length;

  SleqpNewtonData* newton_data;

  SleqpSparseVec* newton_step;

  SleqpSparseVec* cauchy_newton_direction;

  SleqpSparseVec* cauchy_newton_hessian_direction;

  SleqpSparseVec* trial_direction;

  SleqpSparseVec* initial_trial_point;

  SleqpIterate* trial_iterate;

  SleqpAugJacobian* aug_jacobian;

  SleqpDualEstimationData* estimation_data;
  SleqpSparseVec* estimation_residuum;

  SleqpMeritData* merit_data;

  SleqpSparseVec* linear_merit_gradient;

  SleqpSOCData* soc_data;

  SleqpSparseVec* soc_direction;

  SleqpSparseVec* soc_corrected_direction;

  SleqpSparseVec* soc_hessian_direction;

  SleqpSparseVec* initial_soc_trial_point;

  // parameters, adjusted throughout...

  double trust_radius;

  double lp_trust_radius;

  double penalty_parameter;

  // misc

  int iteration;
};


SLEQP_RETCODE sleqp_solver_create(SleqpSolver** star,
                                  SleqpProblem* problem,
                                  SleqpParams* params,
                                  SleqpSparseVec* x)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpSolver* solver = *star;

  const int num_constraints = problem->num_constraints;
  const int num_variables = problem->num_variables;

  solver->problem = problem;
  solver->params = params;
  solver->iteration = 0;

  SLEQP_CALL(sleqp_deriv_checker_create(&solver->deriv_check,
                                        problem,
                                        params));

  const double eps = sleqp_params_get_eps(params);

  assert(sleqp_sparse_vector_valid(x));

  SLEQP_CALL(sleqp_iterate_create(&solver->iterate,
                                  solver->problem,
                                  x));

  SLEQP_CALL(sleqp_sparse_vector_clip(x,
                                      problem->var_lb,
                                      problem->var_ub,
                                      eps,
                                      solver->iterate->x));

  SLEQP_CALL(sleqp_sparse_vector_create(&solver->violation,
                                        num_constraints,
                                        0));

  // TODO: make this generic at a later point...

  int num_lp_variables = num_variables + 2*num_constraints;
  int num_lp_constraints = num_constraints;

  SLEQP_CALL(sleqp_lpi_soplex_create_interface(&solver->lp_interface,
                                               num_lp_variables,
                                               num_lp_constraints));

  SLEQP_CALL(sleqp_cauchy_data_create(&solver->cauchy_data,
                                      problem,
                                      params,
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
                                      problem,
                                      params));

  SLEQP_CALL(sleqp_sparse_vector_create(&solver->newton_step,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&solver->cauchy_newton_direction,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&solver->cauchy_newton_hessian_direction,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&solver->trial_direction,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&solver->initial_trial_point,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_iterate_create(&solver->trial_iterate,
                                  solver->problem,
                                  solver->iterate->x));

  SLEQP_CALL(sleqp_aug_jacobian_create(&solver->aug_jacobian,
                                       problem,
                                       params));

  SLEQP_CALL(sleqp_dual_estimation_data_create(&solver->estimation_data,
                                               problem));

  SLEQP_CALL(sleqp_sparse_vector_create(&solver->estimation_residuum,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_merit_data_create(&solver->merit_data,
                                     problem,
                                     params));

  SLEQP_CALL(sleqp_sparse_vector_create(&solver->linear_merit_gradient,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_soc_data_create(&solver->soc_data,
                                   problem,
                                   params));

  SLEQP_CALL(sleqp_sparse_vector_create(&solver->soc_direction,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&solver->soc_corrected_direction,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&solver->soc_hessian_direction,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&solver->initial_soc_trial_point,
                                        num_variables,
                                        0));


  // initial trust region radii as suggested,
  // penalty parameter as suggested:

  solver->trust_radius = 1.;
  solver->lp_trust_radius = .8 * (solver->trust_radius) * sqrt((double) num_variables);

  solver->penalty_parameter = 10.;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE update_penalty_parameter(SleqpSolver* solver)
{
  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;

  int num_constraints = problem->num_constraints;

  const double violation_tolerance = 1e-8;
  const double min_decrease = .1;
  const int max_increases = 100;

  const double penalty_increase = 10.;

  const double eps = sleqp_params_get_eps(solver->params);

  if(num_constraints == 0)
  {
    return SLEQP_OKAY;
  }

  SleqpCauchyData* cauchy_data = solver->cauchy_data;

  double current_violation;

  sleqp_cauchy_get_violation(cauchy_data, &current_violation);

  current_violation /= num_constraints;

  if(current_violation <= violation_tolerance)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_cauchy_solve(cauchy_data,
                                NULL,
                                solver->penalty_parameter));

  double inf_violation;

  sleqp_cauchy_get_violation(cauchy_data, &inf_violation);

  inf_violation /= num_constraints;

  assert(sleqp_ge(current_violation, inf_violation, eps));

  if(inf_violation <= violation_tolerance)
  {
    for(int i = 0; i < max_increases; ++i)
    {
      solver->penalty_parameter *= penalty_increase;

      SLEQP_CALL(sleqp_cauchy_solve(cauchy_data,
                                    iterate->func_grad,
                                    solver->penalty_parameter));

      double next_violation;

      sleqp_cauchy_get_violation(cauchy_data, &next_violation);

      next_violation /= num_constraints;

      if(next_violation <= violation_tolerance)
      {
        return SLEQP_OKAY;
      }
    }
  }
  else
  {
    if(current_violation - inf_violation <= violation_tolerance)
    {
      // we can't make progress in feasibility, no need for an increase
      return SLEQP_OKAY;
    }

    for(int i = 0; i < max_increases; ++i)
    {
      solver->penalty_parameter *= penalty_increase;

      SLEQP_CALL(sleqp_cauchy_solve(cauchy_data,
                                    iterate->func_grad,
                                    solver->penalty_parameter));

      double next_violation;

      sleqp_cauchy_get_violation(cauchy_data, &next_violation);

      next_violation /= num_constraints;

      if((current_violation - next_violation) >= min_decrease*(current_violation - inf_violation))
      {
        return SLEQP_OKAY;
      }
    }
  }


  return SLEQP_OKAY;
}

static SLEQP_RETCODE update_lp_trust_radius(bool trial_step_accepted,
                                            double trial_direction_infnorm,
                                            double cauchy_step_infnorm,
                                            double cauchy_step_length,
                                            double eps,
                                            double* lp_trust_radius)
{
  if(trial_step_accepted)
  {
    double norm_increase_factor = 1.2;

    trial_direction_infnorm *= norm_increase_factor;
    cauchy_step_infnorm *= norm_increase_factor;

    double scaled_trust_radius = .1 * (*lp_trust_radius);

    double update_lhs = SLEQP_MAX(trial_direction_infnorm,
                                  cauchy_step_infnorm);

    update_lhs = SLEQP_MAX(update_lhs, scaled_trust_radius);

    if(sleqp_eq(cauchy_step_length, 1., eps))
    {
      (*lp_trust_radius) *= 7.;
    }

    *lp_trust_radius = SLEQP_MIN(update_lhs, *lp_trust_radius);

  }
  else
  {
    double half_norm = .5 * trial_direction_infnorm;
    double small_radius = .1 * (*lp_trust_radius);

    double reduced_radius = SLEQP_MAX(half_norm, small_radius);

    *lp_trust_radius = SLEQP_MIN(reduced_radius, *lp_trust_radius);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE update_trust_radius(double reduction_ratio,
                                         bool trial_step_accepted,
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
  else if(trial_step_accepted)
  {
    // stays the same
  }
  else
  {
    *trust_radius = SLEQP_MIN(0.5*(*trust_radius),
                              0.5*direction_norm);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE compute_trial_direction(SleqpSolver* solver,
                                             double* quadratic_value)
{
  const double eps = sleqp_params_get_eps(solver->params);

  SLEQP_CALL(sleqp_merit_linear_gradient(solver->merit_data,
                                         solver->iterate,
                                         solver->cauchy_step,
                                         solver->penalty_parameter,
                                         solver->linear_merit_gradient));

  double prod_chc;
  double prod_che;
  double prod_ehe;

  SLEQP_CALL(sleqp_sparse_vector_dot(solver->cauchy_step,
                                     solver->cauchy_hessian_direction,
                                     &prod_chc));

  SLEQP_CALL(sleqp_sparse_vector_dot(solver->cauchy_step,
                                     solver->cauchy_newton_hessian_direction,
                                     &prod_che));

  SLEQP_CALL(sleqp_sparse_vector_dot(solver->cauchy_newton_direction,
                                     solver->cauchy_newton_hessian_direction,
                                     &prod_ehe));

  // holds the product of the gradient of the quadratic model
  // at alpha = 0 with the step direction
  double gradient_product;

  SLEQP_CALL(sleqp_sparse_vector_dot(solver->cauchy_newton_direction,
                                     solver->linear_merit_gradient,
                                     &gradient_product));

  // holds the value of the quadratic model at alpha = 0
  double zero_func_val = 0.;

  {
    SLEQP_CALL(sleqp_merit_linear(solver->merit_data,
                                  solver->iterate,
                                  solver->cauchy_step,
                                  solver->penalty_parameter,
                                  &zero_func_val));

    zero_func_val += .5 * prod_chc;
  }

  double alpha = 1.;

  // compute the maximum trial step length
  {
    SleqpProblem* problem = solver->problem;

    SLEQP_CALL(sleqp_sparse_vector_add(solver->iterate->x,
                                       solver->cauchy_step,
                                       eps,
                                       solver->trial_direction));

    SLEQP_CALL(sleqp_max_step_length(solver->trial_direction,
                                     solver->cauchy_newton_direction,
                                     problem->var_lb,
                                     problem->var_ub,
                                     &alpha));
  }

  double eta = sleqp_params_get_linesearch_eta(solver->params);
  double tau = sleqp_params_get_linesearch_tau(solver->params);
  double cutoff_threshold = sleqp_params_get_linesearch_cutoff(solver->params);

  int max_it = 10000;
  int it = 0;

  for(it = 0; it < max_it; ++it)
  {
    // holds the value of the quadratic model at alpha
    double alpha_func_val = 0.;

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(solver->cauchy_step,
                                              solver->cauchy_newton_direction,
                                              1.,
                                              alpha,
                                              eps,
                                              solver->trial_direction));

    {
      double temp;

      SLEQP_CALL(sleqp_merit_linear(solver->merit_data,
                                    solver->iterate,
                                    solver->trial_direction,
                                    solver->penalty_parameter,
                                    &temp));

      alpha_func_val += temp;

      alpha_func_val += .5 * prod_chc;

      alpha_func_val += alpha * prod_che;

      alpha_func_val += .5*(alpha*alpha) * prod_ehe;
    }

    (*quadratic_value) = alpha_func_val;

    if(alpha_func_val <= zero_func_val + eta*alpha*gradient_product)
    {
      break;
    }

    // if alpha becomes too small, set it to zero
    if(alpha <= cutoff_threshold)
    {
      alpha = 0.;

      SLEQP_CALL(sleqp_sparse_vector_copy(solver->cauchy_step,
                                          solver->trial_direction));

      break;
    }

    alpha *= tau;
  }

  assert(it != max_it);

  sleqp_log_debug("Cauchy-Newton line search converged after %d iterations (final value: %12e)",
                  it,
                  alpha);

  return SLEQP_OKAY;
}

static bool should_terminate(SleqpSolver* solver)
{
  SleqpIterate* iterate = solver->iterate;
  SleqpProblem* problem = solver->problem;

  const double optimality_residuum = sleqp_sparse_vector_norminf(solver->estimation_residuum);

  double multiplier_norm = 0.;

  {
    multiplier_norm += sleqp_sparse_vector_normsq(iterate->cons_dual);
    multiplier_norm += sleqp_sparse_vector_normsq(iterate->vars_dual);

    multiplier_norm = sqrt(multiplier_norm);
  }

  double slackness_residuum = sleqp_iterate_slackness_residuum(iterate, problem);

  const double optimality_tolerance = sleqp_params_get_optimality_tol(solver->params);

  {
    double residuum = SLEQP_MAX(optimality_residuum, slackness_residuum);

    if(residuum >= optimality_tolerance * (1 + multiplier_norm))
    {
      return false;
    }
  }

  double cons_residuum = sleqp_sparse_vector_norminf(solver->violation);

  return (cons_residuum < optimality_tolerance * (1 + multiplier_norm));
}

static SLEQP_RETCODE compute_trial_point(SleqpSolver* solver,
                                         double* quadratic_value)
{
  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;

  double eps = sleqp_params_get_eps(solver->params);

  double one = 1.;

  // compute Cauchy direction / step and dual estimation
  {

    SLEQP_CALL(sleqp_cauchy_set_iterate(solver->cauchy_data,
                                        iterate,
                                        solver->lp_trust_radius));

    SLEQP_CALL(sleqp_cauchy_solve(solver->cauchy_data,
                                  iterate->func_grad,
                                  solver->penalty_parameter));

    SLEQP_CALL(sleqp_cauchy_get_active_set(solver->cauchy_data,
                                           iterate,
                                           solver->lp_trust_radius));

    SLEQP_CALL(sleqp_aug_jacobian_set_iterate(solver->aug_jacobian,
                                              iterate));

    SLEQP_CALL(sleqp_cauchy_get_direction(solver->cauchy_data,
                                          iterate,
                                          solver->cauchy_direction));

    {
      /*
      bool contained;

      SLEQP_CALL(sleqp_iterate_active_set_contains(iterate,
                                                   problem,
                                                   solver->cauchy_direction,
                                                   eps,
                                                   &contained));

      assert(contained);
      */
    }
    /*
    {
      double product;

      SLEQP_CALL(sleqp_sparse_vector_dot(iterate->func_grad,
                                         solver->cauchy_direction,
                                         &product));

      assert(!sleqp_pos(product, eps));

    }
   */

    SLEQP_CALL(sleqp_dual_estimation_compute(solver->estimation_data,
                                             iterate,
                                             solver->estimation_residuum,
                                             solver->aug_jacobian));

    SLEQP_CALL(sleqp_func_hess_product(problem->func,
                                       &one,
                                       solver->cauchy_direction,
                                       iterate->cons_dual,
                                       solver->cauchy_hessian_direction));

    SLEQP_CALL(sleqp_sparse_vector_copy(solver->cauchy_direction,
                                        solver->cauchy_step));

    SLEQP_CALL(sleqp_cauchy_compute_step(solver->cauchy_data,
                                         iterate,
                                         solver->penalty_parameter,
                                         solver->trust_radius,
                                         solver->cauchy_hessian_direction,
                                         solver->cauchy_step,
                                         &solver->cauchy_step_length));

    SLEQP_CALL(sleqp_sparse_vector_scale(solver->cauchy_hessian_direction,
                                         solver->cauchy_step_length));
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

  {
    double cnorm = sqrt(sleqp_sparse_vector_normsq(solver->cauchy_step));
    double nnorm = sqrt(sleqp_sparse_vector_normsq(solver->newton_step));

    double nprod;

    SLEQP_CALL(sleqp_sparse_vector_dot(solver->iterate->func_grad,
                                       solver->newton_step,
                                       &nprod));

    double objval, lin_term, quad_term;

    SLEQP_CALL(sleqp_merit_linear(solver->merit_data,
                                  solver->iterate,
                                  solver->newton_step,
                                  solver->penalty_parameter,
                                  &lin_term));

    SLEQP_CALL(sleqp_func_hess_product(problem->func,
                                       &one,
                                       solver->newton_step,
                                       iterate->cons_dual,
                                       solver->cauchy_newton_hessian_direction));

    SLEQP_CALL(sleqp_sparse_vector_dot(solver->cauchy_newton_hessian_direction,
                                       solver->newton_step,
                                       &quad_term));

    objval = lin_term + 0.5*quad_term;

    sleqp_log_debug("Cauchy step norm: %e, Newton step norm: %e, " \
                    "Newton grad prod: %e, Newton quad objval: %e ",
                    cnorm,
                    nnorm,
                    nprod,
                    objval);


  }

  // compute Cauchy Newton direction (from Cauchy point to Newton point)
  {
    SLEQP_CALL(sleqp_sparse_vector_add_scaled(solver->newton_step,
                                              solver->cauchy_step,
                                              1.,
                                              -1.,
                                              eps,
                                              solver->cauchy_newton_direction));

    SLEQP_CALL(sleqp_func_hess_product(problem->func,
                                       &one,
                                       solver->cauchy_newton_direction,
                                       iterate->cons_dual,
                                       solver->cauchy_newton_hessian_direction));
  }

  SLEQP_CALL(compute_trial_direction(solver, quadratic_value));

  SLEQP_CALL(sleqp_sparse_vector_add(iterate->x,
                                     solver->trial_direction,
                                     eps,
                                     solver->initial_trial_point));

  SLEQP_CALL(sleqp_sparse_vector_clip(solver->initial_trial_point,
                                      problem->var_lb,
                                      problem->var_ub,
                                      eps,
                                      solver->trial_iterate->x));


  return SLEQP_OKAY;
}

static SLEQP_RETCODE compute_soc_trial_point(SleqpSolver* solver,
                                             double* soc_value)
{
  SleqpProblem* problem = solver->problem;

  SleqpIterate* iterate = solver->iterate;
  SleqpIterate* trial_iterate = solver->trial_iterate;

  SleqpSparseVec* current_point = iterate->x;
  SleqpSparseVec* trial_point = trial_iterate->x;

  const double eps = sleqp_params_get_eps(solver->params);

  double soc_zero_value;

  // the quadratic value at zero is just the linear value
  {
    SLEQP_CALL(sleqp_sparse_vector_clear(solver->soc_direction));

    SLEQP_CALL(sleqp_merit_linear(solver->merit_data,
                                  solver->iterate,
                                  solver->soc_direction,
                                  solver->penalty_parameter,
                                  &soc_zero_value));
  }

  SLEQP_CALL(sleqp_soc_compute(solver->soc_data,
                               solver->aug_jacobian,
                               iterate,
                               trial_iterate,
                               solver->soc_direction));

  double max_step_length = 1.;

  SLEQP_CALL(sleqp_max_step_length(trial_point,
                                   solver->soc_direction,
                                   problem->var_lb,
                                   problem->var_ub,
                                   &max_step_length));

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(solver->trial_direction,
                                            solver->soc_direction,
                                            1.,
                                            max_step_length,
                                            eps,
                                            solver->soc_corrected_direction));

  {
    SLEQP_CALL(sleqp_sparse_vector_add_scaled(current_point,
                                              solver->soc_corrected_direction,
                                              1.,
                                              1.,
                                              eps,
                                              solver->initial_soc_trial_point));

    SLEQP_CALL(sleqp_sparse_vector_clip(solver->initial_soc_trial_point,
                                        problem->var_lb,
                                        problem->var_ub,
                                        eps,
                                        trial_point));
  }

  {
    double soc_lin_value;

    SLEQP_CALL(sleqp_merit_linear(solver->merit_data,
                                  solver->iterate,
                                  solver->soc_corrected_direction,
                                  solver->penalty_parameter,
                                  &soc_lin_value));

    double one = 1.;

    double soc_quad_value;

    SLEQP_CALL(sleqp_func_hess_product(problem->func,
                                       &one,
                                       solver->soc_corrected_direction,
                                       iterate->cons_dual,
                                       solver->soc_hessian_direction));

    SLEQP_CALL(sleqp_sparse_vector_dot(solver->soc_corrected_direction,
                                       solver->soc_hessian_direction,
                                       &soc_quad_value));

    soc_quad_value *= .5;

    *soc_value = soc_lin_value + soc_quad_value;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE set_func_value(SleqpSolver* solver,
                                    SleqpIterate* iterate)
{
  SleqpProblem* problem = solver->problem;

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

  return SLEQP_OKAY;
}

#define HEADER_FORMAT "%10s |%20s |%20s |%20s |%20s |%20s\n"

#define LINE_FORMAT SLEQP_BOLD "%10d " SLEQP_NO_BOLD "|%20e |%20e |%20e |%20e |%20e\n"

static SLEQP_RETCODE print_header()
{
  fprintf(stdout,
          HEADER_FORMAT,
          "iter",
          "funcval",
          "violation",
          "penalty",
          "LP trust radius",
          "EQP trust radius");

  return SLEQP_OKAY;
}

static SLEQP_RETCODE print_line(SleqpSolver* solver)
{
  fprintf(stdout,
          LINE_FORMAT,
          solver->iteration,
          solver->iterate->func_val,
          sleqp_sparse_vector_norminf(solver->violation),
          solver->penalty_parameter,
          solver->lp_trust_radius,
          solver->trust_radius);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE sleqp_perform_iteration(SleqpSolver* solver)
{
  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;
  SleqpIterate* trial_iterate = solver->trial_iterate;

  assert(sleqp_sparse_vector_is_boxed(iterate->x,
                                      problem->var_lb,
                                      problem->var_ub));

  const double eps = sleqp_params_get_eps(solver->params);

  const double accepted_reduction = sleqp_params_get_accepted_reduction(solver->params);

  double quadratic_iterate_value;
  double quadratic_trial_value;

  {
    SLEQP_CALL(sleqp_merit_func(solver->merit_data,
                                iterate,
                                solver->penalty_parameter,
                                &quadratic_iterate_value));
  }

  if(solver->iteration % 25 == 0)
  {
    SLEQP_CALL(print_header());
  }

  SLEQP_CALL(print_line(solver));

  {
    //SLEQP_CALL(sleqp_deriv_check_first_order(solver->deriv_check, iterate));

    //SLEQP_CALL(sleqp_deriv_check_second_order(solver->deriv_check, iterate));
  }

  SLEQP_CALL(compute_trial_point(solver, &quadratic_trial_value));

  double quadratic_reduction = quadratic_iterate_value - quadratic_trial_value;

  assert(!sleqp_neg(quadratic_reduction, eps));

  SLEQP_CALL(set_func_value(solver, trial_iterate));

  SLEQP_CALL(sleqp_func_eval(problem->func,
                             NULL,
                             &(trial_iterate->func_val),
                             NULL,
                             trial_iterate->cons_val,
                             NULL));

  double actual_reduction = 0.;

  double exact_iterate_value = quadratic_iterate_value;

  {
    double exact_trial_value;

    SLEQP_CALL(sleqp_merit_func(solver->merit_data,
                                trial_iterate,
                                solver->penalty_parameter,
                                &exact_trial_value));

    actual_reduction = exact_iterate_value - exact_trial_value;

    sleqp_log_debug("Current merit function value: %e, trial merit function value: %e",
                    exact_iterate_value,
                    exact_trial_value);

  }

  double reduction_ratio = actual_reduction / quadratic_reduction;

  sleqp_log_debug("Reduction ratio: %e, actual: %e, quadratic: %e",
                  reduction_ratio,
                  actual_reduction,
                  quadratic_reduction);

  double direction_norm = sleqp_sparse_vector_normsq(solver->trial_direction);

  direction_norm = sqrt(direction_norm);

  sleqp_log_debug("Trial step norm: %e", direction_norm);

  bool step_accepted = true;

  if(reduction_ratio >= accepted_reduction)
  {
    sleqp_log_debug("Trial step accepted");
  }
  else
  {
    sleqp_log_debug("Trial step rejected");

    step_accepted = false;

    if(problem->num_constraints > 0)
    {
      sleqp_log_debug("Computing second-order correction");

      double soc_quadratic_reduction;

      SLEQP_CALL(set_func_value(solver, iterate));

      SLEQP_CALL(compute_soc_trial_point(solver, &quadratic_trial_value));

      soc_quadratic_reduction = quadratic_iterate_value - quadratic_trial_value;

      // in the SOC case it is not guaranteed that
      // there is a quadratic reduction
      if(sleqp_pos(soc_quadratic_reduction, eps))
      {

        SLEQP_CALL(set_func_value(solver, trial_iterate));

        SLEQP_CALL(sleqp_func_eval(problem->func,
                                   NULL,
                                   &trial_iterate->func_val,
                                   NULL,
                                   trial_iterate->cons_val,
                                   NULL));

        double soc_exact_trial_value;

        SLEQP_CALL(sleqp_merit_func(solver->merit_data,
                                    trial_iterate,
                                    solver->penalty_parameter,
                                    &soc_exact_trial_value));

        double soc_actual_reduction = exact_iterate_value - soc_exact_trial_value;

        double soc_reduction_ratio = soc_actual_reduction / soc_quadratic_reduction;

        sleqp_log_debug("SOC Reduction ratio: %e, actual: %e, quadratic: %e",
                        soc_reduction_ratio,
                        soc_actual_reduction,
                        soc_quadratic_reduction);

        if(soc_reduction_ratio >= accepted_reduction)
        {
          step_accepted = true;
        }
      }

      if(step_accepted)
      {
        sleqp_log_debug("Second-order correction accepted");
      }
      else
      {
        sleqp_log_debug("Second-order correction rejected");
      }
    }
  }

  // update trust radii, penalty parameter
  {
    SLEQP_CALL(update_trust_radius(reduction_ratio,
                                   step_accepted,
                                   direction_norm,
                                   &(solver->trust_radius)));

    double trial_direction_infnorm = sleqp_sparse_vector_norminf(solver->trial_direction);
    double cauchy_step_infnorm = sleqp_sparse_vector_norminf(solver->cauchy_step);

    SLEQP_CALL(update_lp_trust_radius(step_accepted,
                                      trial_direction_infnorm,
                                      cauchy_step_infnorm,
                                      solver->cauchy_step_length,
                                      eps,
                                      &(solver->lp_trust_radius)));

    SLEQP_CALL(update_penalty_parameter(solver));
  }

  // update current iterate

  if(step_accepted)
  {
    // get the remaining data to fill the iterate

    SLEQP_CALL(sleqp_func_eval(problem->func,
                               NULL,
                               NULL,
                               trial_iterate->func_grad,
                               trial_iterate->cons_val,
                               trial_iterate->cons_jac));

    solver->trial_iterate = iterate;
    solver->iterate = trial_iterate;

    SLEQP_CALL(sleqp_get_violation(problem,
                                   iterate,
                                   eps,
                                   solver->violation));
  }
  else
  {
    set_func_value(solver, iterate);
  }

  ++(solver->iteration);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solver_solve(SleqpSolver* solver,
                                 int max_num_iterations)
{
  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;

  const double eps = sleqp_params_get_eps(solver->params);

  solver->iteration = 0;

  sleqp_log_info("Solving a problem with %d variables, %d constraints",
                 problem->num_variables,
                 problem->num_constraints);

  SLEQP_CALL(sleqp_set_and_evaluate(problem, iterate));

  //SLEQP_CALL(sleqp_deriv_check_first_order(solver->deriv_check, iterate));

  SLEQP_CALL(sleqp_get_violation(problem,
                                 iterate,
                                 eps,
                                 solver->violation));

  sleqp_log_info("Initial function value: %f",
                 solver->iterate->func_val);

  solver->status = SLEQP_INVALID;

  for(int i = 0; i < max_num_iterations; ++i)
  {
    SLEQP_CALL(sleqp_perform_iteration(solver));

    if(should_terminate(solver))
    {
      sleqp_log_info("Achieved optimality");
      solver->status = SLEQP_OPTIMAL;
      break;
    }
  }

  if(solver->status != SLEQP_OPTIMAL)
  {
    const double infeasibility = sleqp_sparse_vector_norminf(solver->violation);

    if(sleqp_zero(infeasibility, eps))
    {
      solver->status = SLEQP_FEASIBLE;
    }
    else
    {
      solver->status = SLEQP_INFEASIBLE;
    }

  }

  sleqp_log_info("Primal solution: ");

  SLEQP_CALL(sleqp_sparse_vector_fprintf(solver->iterate->x, stdout));

  sleqp_log_info("Variable multipliers: ");

  SLEQP_CALL(sleqp_sparse_vector_fprintf(solver->iterate->vars_dual, stdout));

  sleqp_log_info("Constraint multipliers: ");

  SLEQP_CALL(sleqp_sparse_vector_fprintf(solver->iterate->cons_dual, stdout));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solver_get_solution(SleqpSolver* solver,
                                        SleqpIterate** iterate)
{
  *iterate = solver->iterate;

  return SLEQP_OKAY;
}

SLEQP_STATUS sleqp_solver_get_status(SleqpSolver* solver)
{
  return solver->status;
}

SLEQP_RETCODE sleqp_solver_free(SleqpSolver** star)
{
  SleqpSolver* solver = *star;

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->initial_soc_trial_point));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->soc_hessian_direction));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->soc_corrected_direction));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->soc_direction));

  SLEQP_CALL(sleqp_soc_data_free(&solver->soc_data));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->linear_merit_gradient));

  SLEQP_CALL(sleqp_merit_data_free(&solver->merit_data));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->estimation_residuum));

  SLEQP_CALL(sleqp_dual_estimation_data_free(&solver->estimation_data));

  SLEQP_CALL(sleqp_aug_jacobian_free(&solver->aug_jacobian));

  SLEQP_CALL(sleqp_iterate_free(&solver->trial_iterate));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->initial_trial_point));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->trial_direction));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->cauchy_newton_hessian_direction));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->cauchy_newton_direction));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->newton_step));

  SLEQP_CALL(sleqp_newton_data_free(&solver->newton_data));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->cauchy_step));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->cauchy_hessian_direction));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->cauchy_direction));

  SLEQP_CALL(sleqp_cauchy_data_free(&solver->cauchy_data));

  SLEQP_CALL(sleqp_lpi_free(&solver->lp_interface));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->violation));

  SLEQP_CALL(sleqp_iterate_free(&solver->iterate));

  SLEQP_CALL(sleqp_deriv_checker_free(&solver->deriv_check));

  sleqp_free(star);

  return SLEQP_OKAY;
}
