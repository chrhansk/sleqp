#include "sleqp_solver.h"

#include <assert.h>
#include <math.h>

#include "sleqp_assert.h"
#include "sleqp_cmp.h"
#include "sleqp_feas.h"
#include "sleqp_mem.h"

#include "sleqp_defs.h"

#include "sleqp_deriv_check.h"

#include "sleqp_aug_jacobian.h"

#include "sleqp_cauchy.h"
#include "sleqp_dual_estimation.h"

#include "sparse/sleqp_sparse_factorization_umfpack.h"

#include "sleqp_newton.h"

#include "sleqp_soc.h"
#include "sleqp_timer.h"

#include "sleqp_iterate.h"
#include "sleqp_linesearch.h"
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

  SleqpSparseVec* cauchy_step;

  SleqpSparseVec* cauchy_hessian_step;

  double cauchy_step_length;

  SleqpSparseVec* multipliers;

  SleqpNewtonData* newton_data;

  SleqpSparseVec* newton_step;

  SleqpSparseVec* newton_hessian_step;

  SleqpSparseVec* trial_step;

  SLEQP_STEPTYPE last_step_type;

  SleqpSparseVec* initial_trial_point;

  SleqpIterate* trial_iterate;

  SleqpSparseFactorization* factorization;

  SleqpAugJacobian* aug_jacobian;

  SleqpDualEstimationData* estimation_data;
  SleqpSparseVec* estimation_residuum;

  SleqpMeritData* merit_data;

  SleqpLineSearchData* linesearch;

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

  if(time_limit != SLEQP_NONE)
  {
    double remaining_time = time_limit - sleqp_timer_elapsed(solver->elapsed_timer);

    remaining_time = SLEQP_MAX(remaining_time, 0.);

    return remaining_time;
  }

  return SLEQP_NONE;
}

static SLEQP_RETCODE set_residuum(SleqpSolver* solver)
{
  const double feas_eps = sleqp_params_get_feasibility_tolerance(solver->params);

  SLEQP_CALL(sleqp_iterate_slackness_residuum(solver->iterate,
                                              solver->problem,
                                              &solver->slackness_residuum));

  SLEQP_CALL(sleqp_iterate_feasibility_residuum(solver->iterate,
                                                solver->problem,
                                                feas_eps,
                                                &solver->feasibility_residuum));

  SLEQP_CALL(sleqp_iterate_stationarity_residuum(solver->iterate,
                                                 solver->problem,
                                                 solver->dense_cache,
                                                 &solver->stationarity_residuum));

  return SLEQP_OKAY;
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

  SleqpProblem* scaled_problem;

  solver->unscaled_problem = problem;

  if(scaling_data)
  {
    solver->scaling_data = scaling_data;

    SLEQP_CALL(sleqp_scaling_capture(solver->scaling_data));

    SLEQP_CALL(sleqp_problem_scaling_create(&solver->problem_scaling,
                                            solver->scaling_data,
                                            problem,
                                            params));

    SLEQP_CALL(sleqp_iterate_create(&solver->unscaled_iterate,
                                    solver->unscaled_problem,
                                    primal));

    SLEQP_CALL(sleqp_iterate_create(&solver->unscaled_trial_iterate,
                                    solver->unscaled_problem,
                                    primal));

    SLEQP_CALL(sleqp_problem_scaling_flush(solver->problem_scaling));

    scaled_problem = sleqp_problem_scaling_get_problem(solver->problem_scaling);
  }
  else
  {
    scaled_problem = problem;
  }

  SleqpFunc* func = scaled_problem->func;

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

    SLEQP_CALL(sleqp_problem_create(&solver->problem,
                                    func,
                                    params,
                                    scaled_problem->var_lb,
                                    scaled_problem->var_ub,
                                    scaled_problem->cons_lb,
                                    scaled_problem->cons_ub));
  }

  SLEQP_CALL(sleqp_timer_create(&solver->elapsed_timer));

  SLEQP_CALL(sleqp_params_capture(params));
  solver->params = params;

  SLEQP_CALL(sleqp_options_capture(options));
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

  int num_lp_variables = num_variables + 2*num_constraints;
  int num_lp_constraints = num_constraints;

  SLEQP_CALL(sleqp_lpi_create_default_interface(&solver->lp_interface,
                                                num_lp_variables,
                                                num_lp_constraints,
                                                params));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->unscaled_violation,
                                              num_constraints));

  SLEQP_CALL(sleqp_cauchy_data_create(&solver->cauchy_data,
                                      solver->problem,
                                      params,
                                      solver->lp_interface));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->cauchy_direction,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->cauchy_step,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->cauchy_hessian_step,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->multipliers,
                                              num_constraints));

  SLEQP_CALL(sleqp_newton_data_create(&solver->newton_data,
                                      solver->problem,
                                      params,
                                      options));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->newton_step,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->newton_hessian_step,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->trial_step,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->initial_trial_point,
                                              num_variables));

  SLEQP_CALL(sleqp_iterate_create(&solver->trial_iterate,
                                  solver->problem,
                                  sleqp_iterate_get_primal(solver->iterate)));

  // create sparse factorization

  SLEQP_CALL(sleqp_sparse_factorization_umfpack_create(&solver->factorization,
                                                       params));

  SLEQP_CALL(sleqp_aug_jacobian_create(&solver->aug_jacobian,
                                       solver->problem,
                                       params,
                                       solver->factorization));

  SLEQP_CALL(sleqp_dual_estimation_data_create(&solver->estimation_data,
                                               solver->problem));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->estimation_residuum,
                                              num_variables));

  SLEQP_CALL(sleqp_merit_data_create(&solver->merit_data,
                                     solver->problem,
                                     params));

  SLEQP_CALL(sleqp_linesearch_create(&solver->linesearch,
                                     solver->problem,
                                     params,
                                     solver->merit_data));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->primal_diff,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->cons_dual_diff,
                                              num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->vars_dual_diff,
                                              num_variables));

  SLEQP_CALL(sleqp_soc_data_create(&solver->soc_data,
                                   solver->problem,
                                   params));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->soc_direction,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->soc_corrected_direction,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->soc_hessian_direction,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->initial_soc_trial_point,
                                              num_variables));


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

  solver->time_limit = SLEQP_NONE;

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

  // sleqp_assert_is_geq(current_violation, inf_violation, eps);

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
                                            double trial_step_infnorm,
                                            double cauchy_step_infnorm,
                                            double cauchy_step_length,
                                            double eps,
                                            double* lp_trust_radius)
{
  if(trial_step_accepted)
  {
    double norm_increase_factor = 1.2;

    trial_step_infnorm *= norm_increase_factor;
    cauchy_step_infnorm *= norm_increase_factor;

    double scaled_trust_radius = .1 * (*lp_trust_radius);

    double update_lhs = SLEQP_MAX(trial_step_infnorm,
                                  cauchy_step_infnorm);

    update_lhs = SLEQP_MAX(update_lhs, scaled_trust_radius);

    if(sleqp_is_eq(cauchy_step_length, 1., eps))
    {
      (*lp_trust_radius) *= 7.;
    }

    *lp_trust_radius = SLEQP_MIN(update_lhs, *lp_trust_radius);

  }
  else
  {
    double half_norm = .5 * trial_step_infnorm;
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

  SLEQP_CALL(sleqp_newton_compute_multipliers(solver->newton_data,
                                              solver->multipliers));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE compute_cauchy_step(SleqpSolver* solver,
                                         double* cauchy_merit_value,
                                         bool quadratic_model,
                                         bool* full_step)
{
  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;

  const double eps = sleqp_params_get_eps(solver->params);
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

    SLEQP_CALL(sleqp_newton_set_iterate(solver->newton_data,
                                        iterate,
                                        solver->aug_jacobian,
                                        solver->trust_radius,
                                        solver->penalty_parameter));

    SLEQP_CALL(estimate_dual_values(solver, iterate));

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

    SLEQP_CALL(sleqp_sparse_vector_copy(solver->cauchy_direction,
                                        solver->cauchy_step));

    if(!quadratic_model)
    {
      (*full_step) = true;

      solver->cauchy_step_length = 1.;

      return SLEQP_OKAY;
    }

    SLEQP_CALL(sleqp_func_hess_prod(problem->func,
                                    &one,
                                    solver->cauchy_direction,
                                    solver->multipliers,
                                    solver->cauchy_hessian_step));

    SLEQP_CALL(sleqp_linesearch_set_iterate(solver->linesearch,
                                            iterate,
                                            solver->penalty_parameter,
                                            solver->trust_radius));

    SLEQP_CALL(sleqp_linesearch_cauchy_step(solver->linesearch,
                                            solver->cauchy_step,
                                            solver->multipliers,
                                            solver->cauchy_hessian_step,
                                            &solver->cauchy_step_length,
                                            cauchy_merit_value));

#if !defined(NDEBUG)

    {
      double actual_quadratic_merit_value;

      double func_dual = 1.;

      SLEQP_CALL(sleqp_merit_quadratic(solver->merit_data,
                                       iterate,
                                       &func_dual,
                                       solver->cauchy_step,
                                       solver->multipliers,
                                       solver->penalty_parameter,
                                       &actual_quadratic_merit_value));

      sleqp_assert_is_eq(*cauchy_merit_value,
                         actual_quadratic_merit_value,
                         eps);
    }

#endif

    (*full_step) = sleqp_is_eq(solver->cauchy_step_length,
                               1.,
                               zero_eps);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE compute_trial_point_simple(SleqpSolver* solver,
                                                double* cauchy_merit_value,
                                                bool quadratic_model,
                                                bool* full_step)
{
  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;

  const double eps = sleqp_params_get_eps(solver->params);
  const double zero_eps = sleqp_params_get_zero_eps(solver->params);

  const double one = 1.;

  SLEQP_CALL(compute_cauchy_step(solver,
                                 cauchy_merit_value,
                                 quadratic_model,
                                 full_step));

  const SleqpSparseVec* trial_step = solver->cauchy_step;

  // Compute quadratic merit value
  {
    SLEQP_CALL(sleqp_merit_linear(solver->merit_data,
                                  solver->iterate,
                                  trial_step,
                                  solver->penalty_parameter,
                                  cauchy_merit_value));
  }

  if(quadratic_model)
  {
    double hessian_prod;

    SLEQP_CALL(sleqp_sparse_vector_dot(trial_step,
                                       solver->cauchy_hessian_step,
                                       &hessian_prod));

    (*cauchy_merit_value) += .5 * hessian_prod;

#if !defined(NDEBUG)

    {
      double actual_quadratic_merit_value;

      double func_dual = 1.;

      SLEQP_CALL(sleqp_merit_quadratic(solver->merit_data,
                                       iterate,
                                       &func_dual,
                                       solver->cauchy_step,
                                       solver->multipliers,
                                       solver->penalty_parameter,
                                       &actual_quadratic_merit_value));

      sleqp_assert_is_eq(*cauchy_merit_value,
                         actual_quadratic_merit_value,
                         eps);
    }

#endif

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
                                                double* trial_merit_value,
                                                bool* full_step)
{
  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;

  const double eps = sleqp_params_get_eps(solver->params);
  const double zero_eps = sleqp_params_get_zero_eps(solver->params);

  const double one = 1.;

  double cauchy_merit_value;

  SLEQP_CALL(compute_cauchy_step(solver,
                                 &cauchy_merit_value,
                                 true,
                                 full_step));

  // compute Newton step
  {
    SLEQP_CALL(sleqp_newton_set_time_limit(solver->newton_data,
                                           remaining_time(solver)));

    SLEQP_CALL(sleqp_newton_compute_step(solver->newton_data,
                                         solver->multipliers,
                                         solver->newton_step));

    SLEQP_CALL(sleqp_func_hess_prod(problem->func,
                                    &one,
                                    solver->newton_step,
                                    solver->multipliers,
                                    solver->newton_hessian_step));
  }

  {
    double step_length;

    SLEQP_CALL(sleqp_linesearch_trial_step(solver->linesearch,
                                           solver->cauchy_step,
                                           solver->cauchy_hessian_step,
                                           cauchy_merit_value,
                                           solver->newton_step,
                                           solver->newton_hessian_step,
                                           solver->multipliers,
                                           solver->trial_step,
                                           &step_length,
                                           trial_merit_value));
  }

#if !defined(NDEBUG)

  {
    double actual_quadratic_merit_value;

    double func_dual = 1.;

    SLEQP_CALL(sleqp_merit_quadratic(solver->merit_data,
                                     iterate,
                                     &func_dual,
                                     solver->trial_step,
                                     solver->multipliers,
                                     solver->penalty_parameter,
                                     &actual_quadratic_merit_value));

    sleqp_assert_is_eq(*trial_merit_value,
                       actual_quadratic_merit_value,
                       eps);
  }

#endif

  SLEQP_CALL(sleqp_sparse_vector_add(sleqp_iterate_get_primal(iterate),
                                     solver->trial_step,
                                     zero_eps,
                                     solver->initial_trial_point));

  SLEQP_CALL(sleqp_sparse_vector_clip(solver->initial_trial_point,
                                      problem->var_lb,
                                      problem->var_ub,
                                      zero_eps,
                                      sleqp_iterate_get_primal(solver->trial_iterate)));

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

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(solver->trial_step,
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
                                    solver->multipliers,
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

#define HEADER_FORMAT "%10s |%14s |%14s |%14s |%14s |%14s |%14s |%14s |%14s |%14s |%14s |%14s |%14s | %18s"

#define LINE_FORMAT SLEQP_FORMAT_BOLD "%10d " SLEQP_FORMAT_RESET "|%14e |%14e |%14e |%14e |%14e |%14s |%14e |%14e |%14e |%14s |%14e |%14e | %18s"

static SLEQP_RETCODE print_header()
{
  sleqp_log_info(HEADER_FORMAT,
                 "Iteration",
                 "Func val",
                 "Feas res",
                 "Slack res",
                 "Stat res",
                 "Penalty",
                 "Working set",
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

  char jac_condition_buf[1024];

  if(aug_jac_condition != SLEQP_NONE)
  {
    snprintf(jac_condition_buf,
             1024,
             "%f",
             aug_jac_condition);
  }
  else
  {
    snprintf(jac_condition_buf,
             1024,
             "-");
  }

  char working_set_buf[1024];

  SleqpWorkingSet* working_set = sleqp_iterate_get_working_set(solver->iterate);

  snprintf(working_set_buf,
           1024,
           "%dv/%dc",
           sleqp_working_set_num_active_vars(working_set),
           sleqp_working_set_num_active_cons(working_set));

  sleqp_log_info(LINE_FORMAT,
                 solver->iteration,
                 sleqp_iterate_get_func_val(solver->unscaled_iterate),
                 solver->feasibility_residuum,
                 solver->slackness_residuum,
                 solver->stationarity_residuum,
                 solver->penalty_parameter,
                 working_set_buf,
                 solver->lp_trust_radius,
                 solver->trust_radius,
                 basis_condition,
                 jac_condition_buf,
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

  solver->dual_diff_norm += sleqp_sparse_vector_norm_sq(solver->cons_dual_diff);
  solver->dual_diff_norm += sleqp_sparse_vector_norm_sq(solver->vars_dual_diff);

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

  const bool quadratic_model = sleqp_options_get_use_quadratic_model(options);

  const bool perform_newton_step = sleqp_options_get_perform_newton_step(options) && quadratic_model;

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

  if(perform_newton_step)
  {
    SLEQP_CALL(compute_trial_point_newton(solver,
                                          &model_trial_value,
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
    const SLEQP_DERIV_CHECK deriv_check = sleqp_options_get_deriv_check(options);

    if(deriv_check & SLEQP_DERIV_CHECK_FIRST)
    {
      SLEQP_CALL(sleqp_deriv_check_first_order(solver->deriv_check, iterate));
    }

    if(deriv_check & SLEQP_DERIV_CHECK_SECOND_EXHAUSTIVE)
    {
      SLEQP_CALL(sleqp_deriv_check_second_order_exhaustive(solver->deriv_check, iterate));
    }
    else if(deriv_check & SLEQP_DERIV_CHECK_SECOND_SIMPLE)
    {
      SLEQP_CALL(sleqp_deriv_check_second_order_simple(solver->deriv_check, iterate));
    }
  }

  {
    SLEQP_CALL(set_residuum(solver));

    // We perform the optimality test wrt. the scaled problem
    if(sleqp_iterate_is_optimal(iterate,
                                solver->params,
                                solver->feasibility_residuum,
                                solver->slackness_residuum,
                                solver->stationarity_residuum))
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

  sleqp_assert_is_geq(model_reduction, 0., zero_eps);

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

  const double trial_step_infnorm = sleqp_sparse_vector_inf_norm(solver->trial_step);
  const double cauchy_step_infnorm = sleqp_sparse_vector_inf_norm(solver->cauchy_step);

  double trial_step_norm = sleqp_sparse_vector_norm(solver->trial_step);

  sleqp_log_debug("Trial step norm: %e", trial_step_norm);

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
      if(sleqp_is_pos(soc_model_reduction, zero_eps))
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
    if(perform_newton_step)
    {
      SLEQP_CALL(update_trust_radius(reduction_ratio,
                                     step_accepted,
                                     trial_step_norm,
                                     &(solver->trust_radius)));
    }

    SLEQP_CALL(update_lp_trust_radius(step_accepted,
                                      trial_step_infnorm,
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
                                      solver->iterate,
                                      solver->trial_iterate,
                                      solver->multipliers));
    }

    if(solver->sr1_data)
    {
      SLEQP_CALL(sleqp_sr1_data_push(solver->sr1_data,
                                     solver->iterate,
                                     solver->trial_iterate,
                                     solver->multipliers));
    }

    // perform simple swaps
    solver->trial_iterate = iterate;
    solver->iterate = trial_iterate;

    SleqpIterate* unscaled_iterate = solver->unscaled_iterate;
    solver->unscaled_iterate = solver->unscaled_trial_iterate;
    solver->unscaled_trial_iterate = unscaled_iterate;

    SLEQP_CALL(sleqp_violation_values(solver->unscaled_problem,
                                      solver->unscaled_iterate,
                                      zero_eps,
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

  SLEQP_CALL(sleqp_violation_values(problem,
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
    if(time_limit != SLEQP_NONE)
    {
      if(solver->elapsed_seconds >= time_limit)
      {
        break;
      }
    }

    if(max_num_iterations != SLEQP_NONE &&
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

  double violation;

  const double feas_eps = sleqp_params_get_feasibility_tolerance(solver->params);

  SLEQP_CALL(sleqp_iterate_feasibility_residuum(solver->unscaled_iterate,
                                                solver->unscaled_problem,
                                                feas_eps,
                                                &violation));

  if(solver->status != SLEQP_OPTIMAL)
  {
    const bool feasible = sleqp_iterate_is_feasible(iterate,
                                                    solver->feasibility_residuum,
                                                    feas_eps);

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
  const double feas_eps = sleqp_params_get_feasibility_tolerance(solver->params);

  SLEQP_CALL(sleqp_iterate_get_violated_constraints(iterate,
                                                    solver->unscaled_problem,
                                                    violated_constraints,
                                                    num_violated_constraints,
                                                    feas_eps));

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

  SLEQP_CALL(sleqp_soc_data_release(&solver->soc_data));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->vars_dual_diff));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->cons_dual_diff));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->primal_diff));

  SLEQP_CALL(sleqp_linesearch_release(&solver->linesearch));

  SLEQP_CALL(sleqp_merit_data_release(&solver->merit_data));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->estimation_residuum));

  SLEQP_CALL(sleqp_dual_estimation_data_free(&solver->estimation_data));

  SLEQP_CALL(sleqp_aug_jacobian_release(&solver->aug_jacobian));

  SLEQP_CALL(sleqp_sparse_factorization_release(&solver->factorization));

  SLEQP_CALL(sleqp_iterate_release(&solver->trial_iterate));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->initial_trial_point));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->trial_step));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->newton_hessian_step));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->newton_step));

  SLEQP_CALL(sleqp_newton_data_release(&solver->newton_data));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->multipliers));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->cauchy_hessian_step));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->cauchy_step));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->cauchy_direction));

  SLEQP_CALL(sleqp_cauchy_data_free(&solver->cauchy_data));

  SLEQP_CALL(sleqp_lpi_free(&solver->lp_interface));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->unscaled_violation));

  SLEQP_CALL(sleqp_iterate_release(&solver->iterate));

  SLEQP_CALL(sleqp_deriv_checker_free(&solver->deriv_check));

  SLEQP_CALL(sleqp_options_release(&solver->options));
  SLEQP_CALL(sleqp_params_release(&solver->params));

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

  SLEQP_CALL(sleqp_problem_free(&solver->problem));

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
