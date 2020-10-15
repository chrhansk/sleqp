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
#include "sleqp_timer.h"

#include "sleqp_iterate.h"
#include "sleqp_merit.h"
#include "sleqp_scale.h"
#include "sleqp_problem_scaling.h"

#include "lp/sleqp_lpi.h"

#include "sleqp_bfgs.h"
#include "sleqp_sr1.h"

struct SleqpSolver
{
  int refcount;

  SleqpProblem* unscaled_problem;

  SleqpScalingData* scaling_data;

  SleqpProblemScaling* problem_scaling;

  SleqpIterate* unscaled_iterate;

  SleqpIterate* unscaled_trial_iterate;

  SleqpProblem* problem;

  SleqpTimer* elapsed_timer;

  SLEQP_STATUS status;

  SleqpParams* params;

  SleqpOptions* options;

  SleqpDerivCheckData* deriv_check;

  SleqpIterate* iterate;

  SleqpSparseVec* unscaled_violation;

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

  SLEQP_STEPTYPE last_step_type;

  SleqpSparseVec* initial_trial_point;

  SleqpIterate* trial_iterate;

  SleqpAugJacobian* aug_jacobian;

  SleqpDualEstimationData* estimation_data;
  SleqpSparseVec* estimation_residuum;

  SleqpMeritData* merit_data;

  SleqpSparseVec* linear_merit_gradient;

  // Primal / dual step lengths

  SleqpSparseVec* primal_diff;

  double primal_diff_norm;

  SleqpSparseVec* cons_dual_diff;

  SleqpSparseVec* vars_dual_diff;

  double dual_diff_norm;

  // SOC related
  SleqpSOCData* soc_data;

  SleqpSparseVec* soc_direction;

  SleqpSparseVec* soc_corrected_direction;

  SleqpSparseVec* soc_hessian_direction;

  SleqpSparseVec* initial_soc_trial_point;

  double* dense_cache;

  // residuum

  double slackness_residuum;

  double stationarity_residuum;

  double feasibility_residuum;

  // BFGS related

  SleqpBFGSData* bfgs_data;

  // SR1 related

  SleqpSR1Data* sr1_data;

  // parameters, adjusted throughout...

  double trust_radius;

  double lp_trust_radius;

  double penalty_parameter;

  // misc

  double elapsed_seconds;

  int iteration;

  double time_limit;

};

static double remaining_time(SleqpSolver* solver)
{
  double time_limit = solver->time_limit;

  if(time_limit != -1)
  {
    double remaining_time = time_limit - sleqp_timer_elapsed(solver->elapsed_timer);

    remaining_time = SLEQP_MAX(remaining_time, 0.);

    return remaining_time;
  }

  return -1;
}

static void set_residuum(SleqpSolver* solver)
{
  solver->slackness_residuum = sleqp_iterate_slackness_residuum(solver->iterate,
                                                                solver->problem);

  solver->feasibility_residuum = sleqp_iterate_feasibility_residuum(solver->iterate,
                                                                    solver->problem);

  solver->stationarity_residuum = sleqp_iterate_stationarity_residuum(solver->iterate,
                                                                      solver->problem,
                                                                      solver->dense_cache);

}


SLEQP_RETCODE sleqp_solver_create(SleqpSolver** star,
                                  SleqpProblem* problem,
                                  SleqpParams* params,
                                  SleqpOptions* options,
                                  SleqpSparseVec* primal,
                                  SleqpScalingData* scaling_data)
{
  assert(sleqp_sparse_vector_valid(primal));

  SLEQP_CALL(sleqp_malloc(star));

  SleqpSolver* solver = *star;

  *solver = (SleqpSolver){0};

  solver->refcount = 1;

  const int num_constraints = problem->num_constraints;
  const int num_variables = problem->num_variables;

  SleqpFunc* func = problem->func;

  {
    const SLEQP_HESSIAN_EVAL hessian_eval = sleqp_options_get_hessian_eval(options);

    const int num_iter = sleqp_options_get_quasi_newton_num_iterates(options);

    if(hessian_eval == SLEQP_HESSIAN_EVAL_SIMPLE_BFGS ||
       hessian_eval == SLEQP_HESSIAN_EVAL_DAMPED_BFGS)
    {
      const bool damped_bfgs = (hessian_eval == SLEQP_HESSIAN_EVAL_DAMPED_BFGS);

      SLEQP_CALL(sleqp_bfgs_data_create(&solver->bfgs_data,
                                        func,
                                        params,
                                        num_iter,
                                        damped_bfgs));

      func = sleqp_bfgs_get_func(solver->bfgs_data);
    }

    if(hessian_eval == SLEQP_HESSIAN_EVAL_SR1)
    {
      SLEQP_CALL(sleqp_sr1_data_create(&solver->sr1_data,
                                       func,
                                       params,
                                       num_iter));

      func = sleqp_sr1_get_func(solver->sr1_data);
    }

    SLEQP_CALL(sleqp_problem_create(&solver->unscaled_problem,
                                    func,
                                    params,
                                    problem->var_lb,
                                    problem->var_ub,
                                    problem->cons_lb,
                                    problem->cons_ub));
  }

  solver->scaling_data = scaling_data;

  if(solver->scaling_data)
  {
    SLEQP_CALL(sleqp_scaling_capture(solver->scaling_data));

    SLEQP_CALL(sleqp_problem_scaling_create(&solver->problem_scaling,
                                            solver->scaling_data,
                                            problem,
                                            solver->params));

    SLEQP_CALL(sleqp_problem_scaling_set_func(solver->problem_scaling,
                                              func));

    SLEQP_CALL(sleqp_iterate_create(&solver->unscaled_iterate,
                                    solver->unscaled_problem,
                                    primal));

    SLEQP_CALL(sleqp_iterate_create(&solver->unscaled_trial_iterate,
                                    solver->unscaled_problem,
                                    primal));

    SLEQP_CALL(sleqp_problem_scaling_flush(solver->problem_scaling));

    solver->problem = sleqp_problem_scaling_get_problem(solver->problem_scaling);
  }
  else
  {
    solver->problem = solver->unscaled_problem;
  }

  SLEQP_CALL(sleqp_timer_create(&solver->elapsed_timer));

  solver->params = params;
  solver->options = options;

  solver->iteration = 0;
  solver->elapsed_seconds = 0.;

  SLEQP_CALL(sleqp_deriv_checker_create(&solver->deriv_check,
                                        solver->problem,
                                        params));

  const double zero_eps = sleqp_params_get_zero_eps(params);

  SLEQP_CALL(sleqp_iterate_create(&solver->iterate,
                                  solver->problem,
                                  primal));

  SLEQP_CALL(sleqp_sparse_vector_clip(primal,
                                      solver->unscaled_problem->var_lb,
                                      solver->unscaled_problem->var_ub,
                                      zero_eps,
                                      sleqp_iterate_get_primal(solver->iterate)));

  if(solver->scaling_data)
  {
    SLEQP_CALL(sleqp_scale_point(solver->scaling_data,
                                 sleqp_iterate_get_primal(solver->iterate)));
  }

  SLEQP_CALL(sleqp_sparse_vector_create(&solver->unscaled_violation,
                                        num_constraints,
                                        0));

  int num_lp_variables = num_variables + 2*num_constraints;
  int num_lp_constraints = num_constraints;

  SLEQP_CALL(sleqp_lpi_create_default_interface(&solver->lp_interface,
                                                num_lp_variables,
                                                num_lp_constraints,
                                                params));

  SLEQP_CALL(sleqp_cauchy_data_create(&solver->cauchy_data,
                                      solver->problem,
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
                                      solver->problem,
                                      params,
                                      options));

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
                                  sleqp_iterate_get_primal(solver->iterate)));

  SLEQP_CALL(sleqp_aug_jacobian_create(&solver->aug_jacobian,
                                       solver->problem,
                                       params));

  SLEQP_CALL(sleqp_dual_estimation_data_create(&solver->estimation_data,
                                               solver->problem));

  SLEQP_CALL(sleqp_sparse_vector_create(&solver->estimation_residuum,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_merit_data_create(&solver->merit_data,
                                     solver->problem,
                                     params));

  SLEQP_CALL(sleqp_sparse_vector_create(&solver->linear_merit_gradient,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&solver->primal_diff,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&solver->cons_dual_diff,
                                        num_constraints,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&solver->vars_dual_diff,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_soc_data_create(&solver->soc_data,
                                   solver->problem,
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


  if(!solver->scaling_data)
  {
    solver->unscaled_iterate = solver->iterate;

    solver->unscaled_trial_iterate = solver->trial_iterate;
  }

  SLEQP_CALL(sleqp_calloc(&solver->dense_cache, SLEQP_MAX(num_variables, num_constraints)));

  // initial trust region radii as suggested,
  // penalty parameter as suggested:

  solver->trust_radius = 1.;
  solver->lp_trust_radius = .8 * (solver->trust_radius) * sqrt((double) num_variables);

  solver->penalty_parameter = 10.;

  solver->time_limit = -1;

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

  sleqp_log_debug("Updating penalty parameter, average violation is %.10e",
                  current_violation);

  if(current_violation <= violation_tolerance)
  {
    sleqp_log_debug("Average violation is already below the tolerance of %.10e",
                    violation_tolerance);

    return SLEQP_OKAY;
  }

  sleqp_log_debug("Resolving linearization to compute minimum average violation");

  SLEQP_CALL(sleqp_cauchy_solve(cauchy_data,
                                NULL,
                                solver->penalty_parameter));

  {
    bool locally_infeasible;

    SLEQP_CALL(sleqp_cauchy_locally_infeasible(cauchy_data,
                                               &locally_infeasible));

    if(locally_infeasible)
    {
      sleqp_log_warn("Current iterate is locally infeasible");
    }
  }

  double inf_violation;

  sleqp_cauchy_get_violation(cauchy_data, &inf_violation);

  inf_violation /= num_constraints;

  sleqp_log_debug("Minimum average violation: %.10e", inf_violation);

  // assert(sleqp_ge(current_violation, inf_violation, eps));

  if(inf_violation <= violation_tolerance)
  {
    sleqp_log_debug("Minimum average violation is below tolerance");

    for(int i = 0; i < max_increases; ++i)
    {
      solver->penalty_parameter *= penalty_increase;

      sleqp_log_debug("Resolving linearization to compute average violation for penalty value %e",
                      solver->penalty_parameter);

      SLEQP_CALL(sleqp_cauchy_solve(cauchy_data,
                                    sleqp_iterate_get_func_grad(iterate),
                                    solver->penalty_parameter));

      double next_violation;

      sleqp_cauchy_get_violation(cauchy_data, &next_violation);

      next_violation /= num_constraints;

      sleqp_log_debug("Average violation for penalty value %e is %.10e",
                      solver->penalty_parameter,
                      next_violation);

      if(next_violation <= violation_tolerance)
      {
        sleqp_log_debug("Average violation is below the tolerance of %e",
                        solver->penalty_parameter);

        return SLEQP_OKAY;
      }
      else
      {
        sleqp_log_debug("Average violation is above the tolerance of %e, continuing",
                        solver->penalty_parameter);
      }
    }
  }
  else
  {
    sleqp_log_debug("Minimum average violation is above tolerance");

    if(current_violation - inf_violation <= violation_tolerance)
    {
      sleqp_log_debug("Cannot make progress towards feasibility, aborting");
      // we can't make progress in feasibility, no need for an increase
      return SLEQP_OKAY;
    }

    for(int i = 0; i < max_increases; ++i)
    {
      solver->penalty_parameter *= penalty_increase;

      sleqp_log_debug("Resolving linearization to compute average violation for penalty value %e",
                      solver->penalty_parameter);

      SLEQP_CALL(sleqp_cauchy_solve(cauchy_data,
                                    sleqp_iterate_get_func_grad(iterate),
                                    solver->penalty_parameter));

      double next_violation;

      sleqp_cauchy_get_violation(cauchy_data, &next_violation);

      next_violation /= num_constraints;

      sleqp_log_debug("Average violation for penalty value %e is %.10e",
                      solver->penalty_parameter,
                      next_violation);

      if((current_violation - next_violation) >= min_decrease*(current_violation - inf_violation))
      {
        sleqp_log_debug("Penalty value of %e achieves sufficiently high reduction in average violation",
                        solver->penalty_parameter);

        return SLEQP_OKAY;
      }
      else
      {
        sleqp_log_debug("Penalty value of %e does not achieve sufficiently high reduction in average violation",
                        solver->penalty_parameter);
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
                                             double* quadratic_model_objective)
{
  const double zero_eps = sleqp_params_get_zero_eps(solver->params);

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

    SLEQP_CALL(sleqp_sparse_vector_add(sleqp_iterate_get_primal(solver->iterate),
                                       solver->cauchy_step,
                                       zero_eps,
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
                                              zero_eps,
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

    (*quadratic_model_objective) = alpha_func_val;

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

static SLEQP_RETCODE estimate_dual_values(SleqpSolver* solver,
                                          SleqpIterate* iterate)
{
  SleqpOptions* options = solver->options;

  SLEQP_DUAL_ESTIMATION_TYPE estimation_type = sleqp_options_get_dual_estimation_type(options);

  if(estimation_type == SLEQP_DUAL_ESTIMATION_TYPE_LSQ)
  {
    SLEQP_CALL(sleqp_dual_estimation_compute(solver->estimation_data,
                                             iterate,
                                             solver->estimation_residuum,
                                             solver->aug_jacobian));
  }
  else
  {
    assert(estimation_type == SLEQP_DUAL_ESTIMATION_TYPE_LP);

    SLEQP_CALL(sleqp_cauchy_get_dual_estimation(solver->cauchy_data,
                                                iterate));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE compute_linear_step(SleqpSolver* solver,
                                         double* model_objective,
                                         bool quadratic_model,
                                         bool* full_step)
{
  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;

  const double zero_eps = sleqp_params_get_zero_eps(solver->params);

  const double one = 1.;

  // compute Cauchy direction / step and dual estimation
  {

    SLEQP_CALL(sleqp_cauchy_set_iterate(solver->cauchy_data,
                                        iterate,
                                        solver->lp_trust_radius));

    SLEQP_CALL(sleqp_cauchy_solve(solver->cauchy_data,
                                  sleqp_iterate_get_func_grad(iterate),
                                  solver->penalty_parameter));

    SLEQP_CALL(sleqp_cauchy_get_working_set(solver->cauchy_data,
                                            iterate));

    SLEQP_CALL(sleqp_aug_jacobian_set_iterate(solver->aug_jacobian,
                                              iterate));

    SLEQP_CALL(sleqp_cauchy_get_direction(solver->cauchy_data,
                                          solver->cauchy_direction));

#if !defined(NDEBUG)

    {
      const double eps = sleqp_params_get_eps(solver->params);

      bool in_working_set = false;

      SLEQP_CALL(sleqp_direction_in_working_set(problem,
                                                iterate,
                                                solver->cauchy_direction,
                                                solver->dense_cache,
                                                eps,
                                                &in_working_set));

      assert(in_working_set);
    }

#endif

    SLEQP_CALL(estimate_dual_values(solver, iterate));

    SLEQP_CALL(sleqp_func_hess_prod(problem->func,
                                    &one,
                                    solver->cauchy_direction,
                                    sleqp_iterate_get_cons_dual(iterate),
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

    (*full_step) = sleqp_eq(solver->cauchy_step_length,
                            1.,
                            zero_eps);

    SLEQP_CALL(sleqp_sparse_vector_scale(solver->cauchy_hessian_direction,
                                         solver->cauchy_step_length));
  }
}

static SLEQP_RETCODE compute_trial_point_simple(SleqpSolver* solver,
                                                double* model_objective,
                                                bool quadratic_model,
                                                bool* full_step)
{
  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;

  const double zero_eps = sleqp_params_get_zero_eps(solver->params);

  const double one = 1.;

  SLEQP_CALL(compute_linear_step(solver,
                                 model_objective,
                                 quadratic_model,
                                 full_step));

  const SleqpSparseVec* trial_direction = solver->cauchy_step;

  // Compute quadratic merit value
  {
    SLEQP_CALL(sleqp_merit_linear(solver->merit_data,
                                  solver->iterate,
                                  trial_direction,
                                  solver->penalty_parameter,
                                  model_objective));
  }

  if(quadratic_model)
  {
    double hessian_prod;

    SLEQP_CALL(sleqp_sparse_vector_dot(trial_direction,
                                       solver->cauchy_hessian_direction,
                                       &hessian_prod));

    (*model_objective) += .5 * hessian_prod;
  }

  SLEQP_CALL(sleqp_sparse_vector_add(sleqp_iterate_get_primal(iterate),
                                     solver->cauchy_step,
                                     zero_eps,
                                     solver->initial_trial_point));

  SLEQP_CALL(sleqp_sparse_vector_clip(solver->initial_trial_point,
                                      problem->var_lb,
                                      problem->var_ub,
                                      zero_eps,
                                      sleqp_iterate_get_primal(solver->trial_iterate)));

}

static SLEQP_RETCODE compute_trial_point_newton(SleqpSolver* solver,
                                                double* model_objective,
                                                bool quadratic_model,
                                                bool* full_step)
{
  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;

  const double zero_eps = sleqp_params_get_zero_eps(solver->params);

  const double one = 1.;

  SLEQP_CALL(compute_linear_step(solver,
                                 model_objective,
                                 quadratic_model,
                                 full_step));

  // compute Newton step
  {
    SLEQP_CALL(sleqp_newton_set_time_limit(solver->newton_data,
                                           remaining_time(solver)));

    SLEQP_CALL(sleqp_newton_compute_step(solver->newton_data,
                                         solver->iterate,
                                         solver->aug_jacobian,
                                         solver->trust_radius,
                                         solver->penalty_parameter,
                                         solver->newton_step));
  }

  {
    double cnorm = sleqp_sparse_vector_norm(solver->cauchy_step);
    double nnorm = sleqp_sparse_vector_norm(solver->newton_step);

    double nprod;

    SLEQP_CALL(sleqp_sparse_vector_dot(sleqp_iterate_get_func_grad(solver->iterate),
                                       solver->newton_step,
                                       &nprod));

    double objval, lin_term, quad_term;

    SLEQP_CALL(sleqp_merit_linear(solver->merit_data,
                                  solver->iterate,
                                  solver->newton_step,
                                  solver->penalty_parameter,
                                  &lin_term));

    SLEQP_CALL(sleqp_func_hess_prod(problem->func,
                                    &one,
                                    solver->newton_step,
                                    sleqp_iterate_get_cons_dual(iterate),
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
                                              zero_eps,
                                              solver->cauchy_newton_direction));

    SLEQP_CALL(sleqp_func_hess_prod(problem->func,
                                    &one,
                                    solver->cauchy_newton_direction,
                                    sleqp_iterate_get_cons_dual(iterate),
                                    solver->cauchy_newton_hessian_direction));
  }

  SLEQP_CALL(compute_trial_direction(solver, model_objective));

  SLEQP_CALL(sleqp_sparse_vector_add(sleqp_iterate_get_primal(iterate),
                                     solver->trial_direction,
                                     zero_eps,
                                     solver->initial_trial_point));

  SLEQP_CALL(sleqp_sparse_vector_clip(solver->initial_trial_point,
                                      problem->var_lb,
                                      problem->var_ub,
                                      zero_eps,
                                      sleqp_iterate_get_primal(solver->trial_iterate)));

  if(!quadratic_model)
  {
    SLEQP_CALL(sleqp_merit_linear(solver->merit_data,
                                  solver->iterate,
                                  solver->trial_direction,
                                  solver->penalty_parameter,
                                  model_objective));
  }


  return SLEQP_OKAY;
}

static SLEQP_RETCODE compute_soc_trial_point(SleqpSolver* solver,
                                             double* soc_value,
                                             bool quadratic_model)
{
  SleqpProblem* problem = solver->problem;

  SleqpIterate* iterate = solver->iterate;
  SleqpIterate* trial_iterate = solver->trial_iterate;

  SleqpSparseVec* current_point = sleqp_iterate_get_primal(iterate);
  SleqpSparseVec* trial_point = sleqp_iterate_get_primal(trial_iterate);

  const double zero_eps = sleqp_params_get_zero_eps(solver->params);

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
                                            zero_eps,
                                            solver->soc_corrected_direction));

  {
    SLEQP_CALL(sleqp_sparse_vector_add_scaled(current_point,
                                              solver->soc_corrected_direction,
                                              1.,
                                              1.,
                                              zero_eps,
                                              solver->initial_soc_trial_point));

    SLEQP_CALL(sleqp_sparse_vector_clip(solver->initial_soc_trial_point,
                                        problem->var_lb,
                                        problem->var_ub,
                                        zero_eps,
                                        trial_point));
  }

  SLEQP_CALL(sleqp_merit_linear(solver->merit_data,
                                solver->iterate,
                                solver->soc_corrected_direction,
                                solver->penalty_parameter,
                                soc_value));

  if(quadratic_model)
  {
    double one = 1.;

    double soc_quad_value;

    SLEQP_CALL(sleqp_func_hess_prod(problem->func,
                                    &one,
                                    solver->soc_corrected_direction,
                                    sleqp_iterate_get_cons_dual(iterate),
                                    solver->soc_hessian_direction));

    SLEQP_CALL(sleqp_sparse_vector_dot(solver->soc_corrected_direction,
                                       solver->soc_hessian_direction,
                                       &soc_quad_value));

    (*soc_value) += .5 * soc_quad_value;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE set_func_value(SleqpSolver* solver,
                                    SleqpIterate* iterate,
                                    SLEQP_VALUE_REASON reason)
{
  SleqpProblem* problem = solver->problem;

  int func_grad_nnz = 0;
  int cons_val_nnz = 0;
  int cons_jac_nnz = 0;

  SLEQP_CALL(sleqp_func_set_value(problem->func,
                                  sleqp_iterate_get_primal(iterate),
                                  reason,
                                  &func_grad_nnz,
                                  &cons_val_nnz,
                                  &cons_jac_nnz));

  SLEQP_CALL(sleqp_sparse_vector_reserve(sleqp_iterate_get_func_grad(iterate),
                                         func_grad_nnz));

  SLEQP_CALL(sleqp_sparse_vector_reserve(sleqp_iterate_get_cons_val(iterate),
                                         cons_val_nnz));

  SLEQP_CALL(sleqp_sparse_matrix_reserve(sleqp_iterate_get_cons_jac(iterate),
                                         cons_jac_nnz));

  return SLEQP_OKAY;
}

#define HEADER_FORMAT "%10s |%14s |%14s |%14s |%14s |%14s |%14s |%14s |%14s |%14s |%14s |%14s | %18s"

#define LINE_FORMAT SLEQP_FORMAT_BOLD "%10d " SLEQP_FORMAT_RESET "|%14e |%14e |%14e |%14e |%14e |%14e |%14e |%14e |%14e |%14e |%14e | %18s"

static SLEQP_RETCODE print_header()
{
  sleqp_log_info(HEADER_FORMAT,
                 "Iteration",
                 "Func val",
                 "Feas res",
                 "Slack res",
                 "Stat res",
                 "Penalty",
                 "LP tr",
                 "EQP tr",
                 "LP cond",
                 "Jac cond",
                 "Primal step",
                 "Dual step",
                 "Step type");

  return SLEQP_OKAY;
}

static SLEQP_RETCODE print_line(SleqpSolver* solver)
{
  bool exact = false;
  double basis_condition, aug_jac_condition;

  SLEQP_CALL(sleqp_lpi_get_basis_condition(solver->lp_interface,
                                           &exact,
                                           &basis_condition));

  SLEQP_CALL(sleqp_aug_jacobian_get_condition_estimate(solver->aug_jacobian,
                                                       &aug_jac_condition));

  const char* steptype_descriptions[] = {
    [SLEQP_STEPTYPE_NONE] = "",
    [SLEQP_STEPTYPE_ACCEPTED] = "Accepted",
    [SLEQP_STEPTYPE_ACCEPTED_FULL] = "Accepted (full)",
    [SLEQP_STEPTYPE_SOC_ACCEPTED] = "SOC accepted",
    [SLEQP_STEPTYPE_REJECTED] = "Rejected"
  };

  sleqp_log_info(LINE_FORMAT,
                 solver->iteration,
                 sleqp_iterate_get_func_val(solver->unscaled_iterate),
                 solver->feasibility_residuum,
                 solver->slackness_residuum,
                 solver->stationarity_residuum,
                 solver->penalty_parameter,
                 solver->lp_trust_radius,
                 solver->trust_radius,
                 basis_condition,
                 aug_jac_condition,
                 solver->primal_diff_norm,
                 solver->dual_diff_norm,
                 steptype_descriptions[solver->last_step_type]);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE compute_step_lengths(SleqpSolver* solver,
                                          SleqpIterate* previous_iterate,
                                          SleqpIterate* iterate)
{
  const double eps = sleqp_params_get_eps(solver->params);

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(sleqp_iterate_get_primal(previous_iterate),
                                            sleqp_iterate_get_primal(iterate),
                                            1.,
                                            -1,
                                            eps,
                                            solver->primal_diff));

  solver->primal_diff_norm = sleqp_sparse_vector_norm(solver->primal_diff);

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(sleqp_iterate_get_cons_dual(previous_iterate),
                                            sleqp_iterate_get_cons_dual(iterate),
                                            1.,
                                            -1,
                                            eps,
                                            solver->cons_dual_diff));

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(sleqp_iterate_get_vars_dual(previous_iterate),
                                            sleqp_iterate_get_vars_dual(iterate),
                                            1.,
                                            -1,
                                            eps,
                                            solver->vars_dual_diff));

  solver->dual_diff_norm = 0.;

  solver->dual_diff_norm += sleqp_sparse_vector_normsq(solver->cons_dual_diff);
  solver->dual_diff_norm += sleqp_sparse_vector_normsq(solver->vars_dual_diff);

  solver->dual_diff_norm = sqrt(solver->dual_diff_norm);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE sleqp_perform_iteration(SleqpSolver* solver,
                                             bool* optimal,
                                             bool has_previous_iterate)
{
  *optimal = false;

  SLEQP_CALL(sleqp_lpi_set_time_limit(solver->lp_interface, remaining_time(solver)));

  const SleqpOptions* options = solver->options;

  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;
  SleqpIterate* trial_iterate = solver->trial_iterate;

  assert(sleqp_sparse_vector_is_boxed(sleqp_iterate_get_primal(iterate),
                                      problem->var_lb,
                                      problem->var_ub));

  const double eps = sleqp_params_get_eps(solver->params);
  const double zero_eps = sleqp_params_get_zero_eps(solver->params);

  if(has_previous_iterate)
  {
    SLEQP_CALL(compute_step_lengths(solver, trial_iterate, iterate));
  }
  else
  {
    solver->primal_diff_norm = 0.;
    solver->dual_diff_norm = 0.;
  }

  const double accepted_reduction = sleqp_params_get_accepted_reduction(solver->params);

  {
    const SLEQP_DERIV_CHECK deriv_check = sleqp_options_get_deriv_check(options);

    if(deriv_check & SLEQP_DERIV_CHECK_FIRST)
    {
      SLEQP_CALL(sleqp_deriv_check_first_order(solver->deriv_check, iterate));
    }

    if(deriv_check & SLEQP_DERIV_CHECK_SEC)
    {
      SLEQP_CALL(sleqp_deriv_check_second_order(solver->deriv_check, iterate));
    }
  }

  double exact_iterate_value, model_iterate_value;

  {
    SLEQP_CALL(sleqp_merit_func(solver->merit_data,
                                iterate,
                                solver->penalty_parameter,
                                &exact_iterate_value));

    model_iterate_value = exact_iterate_value;
  }

  double model_trial_value;

  bool full_step;

  const bool quadratic_model = sleqp_options_get_use_quadratic_model(options);

  if(sleqp_options_get_perform_newton_step(options))
  {
    SLEQP_CALL(compute_trial_point_newton(solver,
                                          &model_trial_value,
                                          quadratic_model,
                                          &full_step));
  }
  else
  {
    SLEQP_CALL(compute_trial_point_simple(solver,
                                          &model_trial_value,
                                          quadratic_model,
                                          &full_step));
  }



  {
    set_residuum(solver);

    const double optimality_tolerance = sleqp_params_get_optimality_tolerance(solver->params);

    // We perform the optimality test wrt. the scaled problem
    if(sleqp_iterate_is_optimal(iterate,
                                solver->feasibility_residuum,
                                solver->slackness_residuum,
                                solver->stationarity_residuum,
                                optimality_tolerance))
    {
      *optimal = true;
    }
  }

  if(solver->iteration % 25 == 0)
  {
    SLEQP_CALL(print_header());
  }

  SLEQP_CALL(print_line(solver));

  if(*optimal)
  {
    return SLEQP_OKAY;
  }

  double model_reduction = model_iterate_value - model_trial_value;

  assert(!sleqp_neg(model_reduction, zero_eps));

  SLEQP_CALL(set_func_value(solver, trial_iterate, SLEQP_VALUE_REASON_TRYING_ITERATE));

  double func_val;

  SLEQP_CALL(sleqp_func_eval(problem->func,
                             NULL,
                             &func_val,
                             NULL,
                             sleqp_iterate_get_cons_val(trial_iterate),
                             NULL));

  sleqp_iterate_set_func_val(trial_iterate, func_val);

  double actual_reduction = 0.;

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

  double reduction_ratio = 1.;

  if(actual_reduction != model_reduction)
  {
    reduction_ratio = actual_reduction / model_reduction;
  }

  sleqp_log_debug("Reduction ratio: %e, actual: %e, predicted: %e",
                  reduction_ratio,
                  actual_reduction,
                  model_reduction);

  const double trial_direction_infnorm = sleqp_sparse_vector_norminf(solver->trial_direction);
  const double cauchy_step_infnorm = sleqp_sparse_vector_norminf(solver->cauchy_step);

  double trial_direction_norm = sleqp_sparse_vector_norm(solver->trial_direction);

  sleqp_log_debug("Trial step norm: %e", trial_direction_norm);

  bool step_accepted = true;
  bool soc_step_accepted = false;

  solver->last_step_type = SLEQP_STEPTYPE_REJECTED;

  if(reduction_ratio >= accepted_reduction)
  {
    sleqp_log_debug("Trial step accepted");

    if(full_step)
    {
      solver->last_step_type = SLEQP_STEPTYPE_ACCEPTED_FULL;
    }
    else
    {
      solver->last_step_type = SLEQP_STEPTYPE_ACCEPTED;
    }
  }
  else
  {
    sleqp_log_debug("Trial step rejected");

    step_accepted = false;

    if((problem->num_constraints > 0) && sleqp_options_get_perform_soc(options))
    {
      sleqp_log_debug("Computing second-order correction");

      double soc_model_reduction;

      // SLEQP_CALL(set_func_value(solver, iterate));

      SLEQP_CALL(compute_soc_trial_point(solver, &model_trial_value, quadratic_model));

      soc_model_reduction = model_iterate_value - model_trial_value;

      // in the SOC case it is not guaranteed that
      // there is a quadratic reduction
      if(sleqp_pos(soc_model_reduction, zero_eps))
      {

        SLEQP_CALL(set_func_value(solver,
                                  trial_iterate,
                                  SLEQP_VALUE_REASON_TRYING_SOC_ITERATE));

        double func_val;

        SLEQP_CALL(sleqp_func_eval(problem->func,
                                   NULL,
                                   &func_val,
                                   NULL,
                                   sleqp_iterate_get_cons_val(trial_iterate),
                                   NULL));

        SLEQP_CALL(sleqp_iterate_set_func_val(trial_iterate, func_val));

        double soc_exact_trial_value;

        SLEQP_CALL(sleqp_merit_func(solver->merit_data,
                                    trial_iterate,
                                    solver->penalty_parameter,
                                    &soc_exact_trial_value));

        double soc_actual_reduction = exact_iterate_value - soc_exact_trial_value;

        double soc_reduction_ratio = soc_actual_reduction / soc_model_reduction;

        sleqp_log_debug("SOC Reduction ratio: %e, actual: %e, predicted: %e",
                        soc_reduction_ratio,
                        soc_actual_reduction,
                        soc_model_reduction);

        if(soc_reduction_ratio >= accepted_reduction)
        {
          soc_step_accepted = true;
        }
      }

      if(soc_step_accepted)
      {
        solver->last_step_type = SLEQP_STEPTYPE_SOC_ACCEPTED;
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
    if(sleqp_options_get_perform_newton_step(options))
    {
      SLEQP_CALL(update_trust_radius(reduction_ratio,
                                     step_accepted,
                                     trial_direction_norm,
                                     &(solver->trust_radius)));
    }

    SLEQP_CALL(update_lp_trust_radius(step_accepted,
                                      trial_direction_infnorm,
                                      cauchy_step_infnorm,
                                      solver->cauchy_step_length,
                                      zero_eps,
                                      &(solver->lp_trust_radius)));

    SLEQP_CALL(update_penalty_parameter(solver));
  }

  // update current iterate

  if(step_accepted || soc_step_accepted)
  {
    SLEQP_CALL(set_func_value(solver,
                              trial_iterate,
                              SLEQP_VALUE_REASON_ACCEPTED_ITERATE));

    // get the remaining data to fill the iterate

    SLEQP_CALL(sleqp_func_eval(problem->func,
                               NULL,
                               NULL,
                               sleqp_iterate_get_func_grad(trial_iterate),
                               sleqp_iterate_get_cons_val(trial_iterate),
                               sleqp_iterate_get_cons_jac(trial_iterate)));

    // ensure that the unscaled iterate is kept up to date
    if(solver->scaling_data)
    {
      SLEQP_CALL(sleqp_iterate_copy(trial_iterate,
                                    solver->unscaled_trial_iterate));

      SLEQP_CALL(sleqp_unscale_iterate(solver->scaling_data,
                                       solver->unscaled_trial_iterate));
    }

    if(solver->bfgs_data)
    {

      SLEQP_CALL(sleqp_bfgs_data_push(solver->bfgs_data,
                                      solver->unscaled_iterate,
                                      solver->unscaled_trial_iterate));
    }

    if(solver->sr1_data)
    {
      SLEQP_CALL(sleqp_sr1_data_push(solver->sr1_data,
                                     solver->unscaled_iterate,
                                     solver->unscaled_trial_iterate));
    }

    // perform simple swaps
    solver->trial_iterate = iterate;
    solver->iterate = trial_iterate;

    SleqpIterate* unscaled_iterate = solver->unscaled_iterate;
    solver->unscaled_iterate = solver->unscaled_trial_iterate;
    solver->unscaled_trial_iterate = unscaled_iterate;

    SLEQP_CALL(sleqp_get_violation(solver->unscaled_problem,
                                   solver->unscaled_iterate,
                                   eps,
                                   solver->unscaled_violation));
  }
  else
  {
    set_func_value(solver, iterate, SLEQP_VALUE_REASON_REJECTED_ITERATE);
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solver_solve(SleqpSolver* solver,
                                 int max_num_iterations,
                                 double time_limit)
{
  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;

  const double eps = sleqp_params_get_eps(solver->params);
  const double zero_eps = sleqp_params_get_zero_eps(solver->params);

  sleqp_log_info("Solving a problem with %d variables, %d constraints",
                 problem->num_variables,
                 problem->num_constraints);

  SLEQP_CALL(sleqp_set_and_evaluate(problem,
                                    iterate,
                                    SLEQP_VALUE_REASON_INIT));

  // ensure that the unscaled iterate is initialized
  if(solver->scaling_data)
  {
    SLEQP_CALL(sleqp_iterate_copy(iterate,
                                  solver->unscaled_iterate));

    SLEQP_CALL(sleqp_unscale_iterate(solver->scaling_data,
                                     solver->unscaled_iterate));
  }

  SLEQP_CALL(sleqp_get_violation(problem,
                                 iterate,
                                 eps,
                                 solver->unscaled_violation));

  solver->status = SLEQP_INVALID;

  SLEQP_CALL(sleqp_timer_reset(solver->elapsed_timer));

  solver->time_limit = time_limit;

  solver->iteration = 0;
  solver->elapsed_seconds = 0.;
  solver->last_step_type = SLEQP_STEPTYPE_NONE;

  const double deadpoint_bound = sleqp_params_get_deadpoint_bound(solver->params);
  bool reached_deadpoint = false;
  bool has_previous_iterate = false;

  // main solving loop
  while(true)
  {
    if(time_limit != -1)
    {
      if(solver->elapsed_seconds >= time_limit)
      {
        break;
      }
    }

    if(max_num_iterations != -1 &&
       solver->iteration >= max_num_iterations)
    {
      break;
    }

    bool optimal;

    SLEQP_CALL(sleqp_timer_start(solver->elapsed_timer));

    SLEQP_CALL(sleqp_perform_iteration(solver, &optimal, has_previous_iterate));

    switch(solver->last_step_type)
    {
    case SLEQP_STEPTYPE_ACCEPTED:
    case SLEQP_STEPTYPE_ACCEPTED_FULL:
    case SLEQP_STEPTYPE_SOC_ACCEPTED:
      has_previous_iterate = true;
      break;
    default:
      break;
    }

    ++solver->iteration;

    SLEQP_CALL(sleqp_timer_stop(solver->elapsed_timer));

    solver->elapsed_seconds = sleqp_timer_get_ttl(solver->elapsed_timer);

    if(solver->lp_trust_radius <= deadpoint_bound ||
       solver->trust_radius <= deadpoint_bound)
    {
      reached_deadpoint = true;
      break;
    }

    if(optimal)
    {
      sleqp_log_debug("Achieved optimality");
      solver->status = SLEQP_OPTIMAL;
      break;
    }
  }

  if(reached_deadpoint)
  {
    sleqp_log_warn("Reached dead point");
  }

  const double violation = sleqp_iterate_feasibility_residuum(solver->unscaled_iterate,
                                                              solver->unscaled_problem);

  if(solver->status != SLEQP_OPTIMAL)
  {
    const double tolerance = sleqp_params_get_optimality_tolerance(solver->params);

    const bool feasible = sleqp_iterate_is_feasible(iterate,
                                                    solver->feasibility_residuum,
                                                    tolerance);

    if(feasible)
    {
      solver->status = SLEQP_FEASIBLE;
    }
    else
    {
      solver->status = SLEQP_INFEASIBLE;
    }

  }

  const char* descriptions[] = {
    [SLEQP_FEASIBLE] = SLEQP_FORMAT_BOLD SLEQP_FORMAT_YELLOW "feasible" SLEQP_FORMAT_RESET,
    [SLEQP_OPTIMAL] = SLEQP_FORMAT_BOLD SLEQP_FORMAT_GREEN "optimal" SLEQP_FORMAT_RESET,
    [SLEQP_INFEASIBLE] = SLEQP_FORMAT_BOLD SLEQP_FORMAT_RED "infeasible" SLEQP_FORMAT_RESET,
    [SLEQP_INVALID] = SLEQP_FORMAT_BOLD SLEQP_FORMAT_RED "Invalid" SLEQP_FORMAT_RESET
  };

  SleqpFunc* func = solver->unscaled_problem->func;

  sleqp_log_info(SLEQP_FORMAT_BOLD "       Solution status: %s" SLEQP_FORMAT_RESET,
                 descriptions[solver->status]);

  sleqp_log_info(SLEQP_FORMAT_BOLD "       Objective value: %e" SLEQP_FORMAT_RESET,
                 sleqp_iterate_get_func_val(solver->unscaled_iterate));

  sleqp_log_info(SLEQP_FORMAT_BOLD "             Violation: %e" SLEQP_FORMAT_RESET,
                 violation);

  sleqp_log_info(                  "            Iterations: %d", solver->iteration);

  SleqpTimer* func_timer = sleqp_func_get_eval_timer(func);

  SleqpTimer* hess_timer = sleqp_func_get_hess_timer(func);

  SleqpTimer* lp_timer = sleqp_lpi_get_solve_timer(solver->lp_interface);

  sleqp_log_info("  Function evaluations: %4d (%fs avg)",
                 sleqp_timer_get_num_runs(func_timer),
                 sleqp_timer_get_avg(func_timer));

  sleqp_log_info("   Hessian evaluations: %4d (%fs avg)",
                 sleqp_timer_get_num_runs(hess_timer),
                 sleqp_timer_get_avg(hess_timer));

  sleqp_log_info("         LP iterations: %4d (%fs avg)",
                 sleqp_timer_get_num_runs(lp_timer),
                 sleqp_timer_get_avg(lp_timer));

  sleqp_log_info("          Solving time: %.2fs", solver->elapsed_seconds);

  if(solver->status == SLEQP_INFEASIBLE)
  {
    sleqp_log_info("Violations: ");

    for(int index = 0; index < solver->unscaled_violation->nnz; ++index)
    {
      sleqp_log_info("(%d) = %e",
                     solver->unscaled_violation->indices[index],
                     solver->unscaled_violation->data[index]);
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solver_get_solution(SleqpSolver* solver,
                                        SleqpIterate** iterate)
{
  (*iterate) = solver->unscaled_iterate;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solver_get_violated_constraints(SleqpSolver* solver,
                                                    SleqpIterate* iterate,
                                                    int* violated_constraints,
                                                    int* num_violated_constraints)
{
  const double tolerance = sleqp_params_get_optimality_tolerance(solver->params);

  SLEQP_CALL(sleqp_iterate_get_violated_constraints(iterate,
                                                    solver->unscaled_problem,
                                                    tolerance,
                                                    violated_constraints,
                                                    num_violated_constraints));

  return SLEQP_OKAY;
}

SLEQP_STATUS sleqp_solver_get_status(SleqpSolver* solver)
{
  return solver->status;
}

int sleqp_solver_get_iterations(SleqpSolver* solver)
{
  return solver->iteration;
}

double sleqp_solver_get_elapsed_seconds(SleqpSolver* solver)
{
  return solver->elapsed_seconds;
}

static SLEQP_RETCODE solver_free(SleqpSolver** star)
{
  SleqpSolver* solver = *star;

  if(!solver)
  {
    return SLEQP_OKAY;
  }

  sleqp_free(&solver->dense_cache);

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->initial_soc_trial_point));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->soc_hessian_direction));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->soc_corrected_direction));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->soc_direction));

  SLEQP_CALL(sleqp_soc_data_free(&solver->soc_data));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->linear_merit_gradient));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->vars_dual_diff));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->cons_dual_diff));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->primal_diff));

  SLEQP_CALL(sleqp_merit_data_free(&solver->merit_data));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->estimation_residuum));

  SLEQP_CALL(sleqp_dual_estimation_data_free(&solver->estimation_data));

  SLEQP_CALL(sleqp_aug_jacobian_free(&solver->aug_jacobian));

  SLEQP_CALL(sleqp_iterate_release(&solver->trial_iterate));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->initial_trial_point));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->trial_direction));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->cauchy_newton_hessian_direction));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->cauchy_newton_direction));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->newton_step));

  SLEQP_CALL(sleqp_newton_data_release(&solver->newton_data));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->cauchy_step));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->cauchy_hessian_direction));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->cauchy_direction));

  SLEQP_CALL(sleqp_cauchy_data_free(&solver->cauchy_data));

  SLEQP_CALL(sleqp_lpi_free(&solver->lp_interface));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->unscaled_violation));

  SLEQP_CALL(sleqp_iterate_release(&solver->iterate));

  SLEQP_CALL(sleqp_deriv_checker_free(&solver->deriv_check));

  SLEQP_CALL(sleqp_timer_free(&solver->elapsed_timer));

  if(solver->scaling_data)
  {
    SLEQP_CALL(sleqp_iterate_release(&solver->unscaled_trial_iterate));

    SLEQP_CALL(sleqp_iterate_release(&solver->unscaled_iterate));
  }
  else
  {
    solver->unscaled_iterate = NULL;
  }

  SLEQP_CALL(sleqp_problem_scaling_release(&solver->problem_scaling));

  SLEQP_CALL(sleqp_scaling_release(&solver->scaling_data));

  SLEQP_CALL(sleqp_problem_free(&solver->unscaled_problem));

  SLEQP_CALL(sleqp_sr1_data_release(&solver->sr1_data));

  SLEQP_CALL(sleqp_bfgs_data_release(&solver->bfgs_data));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solver_capture(SleqpSolver* solver)
{
  ++solver->refcount;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solver_release(SleqpSolver** star)
{
  SleqpSolver* solver = *star;

  if(!solver)
  {
    return SLEQP_OKAY;
  }

  if(--solver->refcount == 0)
  {
    SLEQP_CALL(solver_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
