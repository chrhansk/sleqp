#include "solver.h"

#include <math.h>

#include "cmp.h"
#include "fail.h"
#include "penalty.h"

static SLEQP_RETCODE set_residuum(SleqpSolver* solver)
{
  SLEQP_CALL(sleqp_iterate_slackness_residuum(solver->problem,
                                              solver->iterate,
                                              &solver->slackness_residuum));

  SLEQP_CALL(sleqp_iterate_feasibility_residuum(solver->problem,
                                                solver->iterate,
                                                &solver->feasibility_residuum));

  SLEQP_CALL(sleqp_iterate_stationarity_residuum(solver->problem,
                                                 solver->iterate,
                                                 solver->dense_cache,
                                                 &solver->stationarity_residuum));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE check_derivative(SleqpSolver* solver)
{
  SleqpOptions* options = solver->options;
  SleqpIterate* iterate = solver->iterate;

  const SLEQP_DERIV_CHECK deriv_check = sleqp_options_get_int(options,
                                                              SLEQP_OPTION_INT_DERIV_CHECK);

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

  return SLEQP_OKAY;
}

static SLEQP_RETCODE compute_step_lengths(SleqpSolver* solver)
{
  SleqpIterate* iterate = solver->iterate;
  SleqpIterate* trial_iterate = solver->trial_iterate;

  const double zero_eps = sleqp_params_get(solver->params, SLEQP_PARAM_ZERO_EPS);

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(sleqp_iterate_get_primal(iterate),
                                            sleqp_iterate_get_primal(trial_iterate),
                                            1.,
                                            -1,
                                            zero_eps,
                                            solver->primal_diff));

  solver->primal_diff_norm = sleqp_sparse_vector_norm(solver->primal_diff);

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(sleqp_iterate_get_cons_dual(iterate),
                                            sleqp_iterate_get_cons_dual(trial_iterate),
                                            1.,
                                            -1,
                                            zero_eps,
                                            solver->cons_dual_diff));

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(sleqp_iterate_get_vars_dual(iterate),
                                            sleqp_iterate_get_vars_dual(trial_iterate),
                                            1.,
                                            -1,
                                            zero_eps,
                                            solver->vars_dual_diff));

  solver->dual_diff_norm = 0.;

  solver->dual_diff_norm += sleqp_sparse_vector_norm_sq(solver->cons_dual_diff);
  solver->dual_diff_norm += sleqp_sparse_vector_norm_sq(solver->vars_dual_diff);

  solver->dual_diff_norm = sqrt(solver->dual_diff_norm);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solver_perform_iteration(SleqpSolver* solver,
                                             bool* optimal)
{
  *optimal = false;

  SLEQP_CALL(sleqp_lpi_set_time_limit(solver->lp_interface,
                                      sleqp_solver_remaining_time(solver)));

  const SleqpOptions* options = solver->options;

  const bool quadratic_model = sleqp_options_get_bool(options,
                                                      SLEQP_OPTION_BOOL_USE_QUADRATIC_MODEL);

  const bool perform_newton_step = quadratic_model &&
    sleqp_options_get_bool(options, SLEQP_OPTION_BOOL_PERFORM_NEWTON_STEP);

  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;
  SleqpIterate* trial_iterate = solver->trial_iterate;

  const int num_constraints = sleqp_problem_num_constraints(problem);

  assert(sleqp_sparse_vector_is_boxed(sleqp_iterate_get_primal(iterate),
                                      sleqp_problem_var_lb(problem),
                                      sleqp_problem_var_ub(problem)));

  const double eps = sleqp_params_get(solver->params,
                                      SLEQP_PARAM_EPS);

  const double zero_eps = sleqp_params_get(solver->params,
                                           SLEQP_PARAM_ZERO_EPS);

  const double accepted_reduction = sleqp_params_get(solver->params,
                                                     SLEQP_PARAM_ACCEPTED_REDUCTION);

  double exact_iterate_value, model_iterate_value;

  {
    SLEQP_CALL(sleqp_merit_func(solver->merit_data,
                                iterate,
                                solver->penalty_parameter,
                                &exact_iterate_value));

    solver->current_merit_value = model_iterate_value = exact_iterate_value;
  }

  SLEQP_CALL(set_residuum(solver));

  if(solver->iteration == 0)
  {
    SLEQP_CALL(sleqp_solver_print_initial_line(solver));
  }

  double model_trial_value;

  bool full_step;

  // Derivative check
  SLEQP_CALL(check_derivative(solver));

  // Optimality check with respect to scaled problem
  {
    if(sleqp_iterate_is_optimal(iterate,
                                solver->params,
                                solver->feasibility_residuum,
                                solver->slackness_residuum,
                                solver->stationarity_residuum))
    {
      *optimal = true;
    }
  }

  if(*optimal)
  {
    return SLEQP_OKAY;
  }

  // Step computation
  if(perform_newton_step)
  {
    SLEQP_CALL(sleqp_solver_compute_trial_point_newton(solver,
                                                       &model_trial_value,
                                                       &full_step));
  }
  else
  {
    SLEQP_CALL(sleqp_solver_compute_trial_point_simple(solver,
                                                       &model_trial_value,
                                                       quadratic_model,
                                                       &full_step));
  }

  SLEQP_CALL(compute_step_lengths(solver));

  double model_reduction = model_iterate_value - model_trial_value;

  sleqp_assert_is_geq(model_reduction, 0., zero_eps);

  SLEQP_CALL(sleqp_solver_set_func_value(solver,
                                         trial_iterate,
                                         SLEQP_VALUE_REASON_TRYING_ITERATE));

  double func_val;

  SLEQP_CALL(sleqp_problem_eval(problem,
                                NULL,
                                &func_val,
                                NULL,
                                sleqp_iterate_get_cons_val(trial_iterate),
                                NULL));

  SLEQP_CALL(sleqp_iterate_set_func_val(trial_iterate, func_val));

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

  solver->boundary_step = sleqp_is_geq(trial_step_norm,
                                       solver->trust_radius,
                                       eps);

  bool step_accepted = true;

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

    const bool perform_soc = sleqp_options_get_bool(options,
                                                    SLEQP_OPTION_BOOL_PERFORM_SOC);

    if((num_constraints > 0) && perform_soc)
    {
      sleqp_log_debug("Computing second-order correction");

      SLEQP_CALL(sleqp_solver_compute_trial_point_soc(solver));

      solver->boundary_step = sleqp_is_geq(solver->soc_step_norm,
                                           solver->trust_radius,
                                           eps);

      SLEQP_CALL(sleqp_solver_set_func_value(solver,
                                             trial_iterate,
                                             SLEQP_VALUE_REASON_TRYING_SOC_ITERATE));

      double func_val;

      SLEQP_CALL(sleqp_problem_eval(problem,
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

      const double soc_actual_reduction = exact_iterate_value - soc_exact_trial_value;

      reduction_ratio = 1.;

      // The denominator of the SOC reduction ratio is the quadratic reduction
      // with respect to the original (not the corrected) trial step
      if(soc_actual_reduction != model_reduction)
      {
        reduction_ratio = soc_actual_reduction / model_reduction;
      }

      sleqp_log_debug("SOC Reduction ratio: %e, actual: %e, predicted: %e",
                      reduction_ratio,
                      soc_actual_reduction,
                      model_reduction);

      if(reduction_ratio >= accepted_reduction)
      {
        step_accepted = true;
      }

      if(step_accepted)
      {
        solver->last_step_type = SLEQP_STEPTYPE_ACCEPTED_SOC;
        sleqp_log_debug("Second-order correction accepted");
      }
      else
      {
        sleqp_log_debug("Second-order correction rejected");
      }
    }
  }

  ++solver->iteration;

  if(solver->iteration % 25 == 0)
  {
    SLEQP_CALL(sleqp_solver_print_header(solver));
  }

  SLEQP_CALL(sleqp_solver_print_line(solver));

  // update trust radii, penalty parameter
  {
    if(perform_newton_step)
    {
      SLEQP_CALL(sleqp_solver_update_trust_radius(solver,
                                                  reduction_ratio,
                                                  step_accepted,
                                                  trial_step_norm));
    }

    SLEQP_CALL(sleqp_solver_update_lp_trust_radius(solver,
                                                   step_accepted,
                                                   trial_step_infnorm,
                                                   cauchy_step_infnorm,
                                                   solver->cauchy_step_length,
                                                   zero_eps,
                                                   &(solver->lp_trust_radius)));

    SLEQP_CALL(sleqp_update_penalty(problem,
                                    iterate,
                                    solver->cauchy_data,
                                    &(solver->penalty_parameter)));
  }

  // update current iterate

  if(step_accepted)
  {
    SLEQP_CALL(sleqp_solver_accept_step(solver));
  }
  else
  {
    SLEQP_CALL(sleqp_solver_reject_step(solver));
  }

  SLEQP_CALLBACK_HANDLER_EXECUTE(solver->callback_handlers[SLEQP_SOLVER_EVENT_PERFORMED_ITERATION],
                                 SLEQP_PERFORMED_ITERATION,
                                 solver);

  return SLEQP_OKAY;
}