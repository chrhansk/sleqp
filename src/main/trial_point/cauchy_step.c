#include "trial_point.h"

#include "cmp.h"
#include "fail.h"
#include "penalty.h"

const double allowed_dual_factor = 1000.;
const double allowed_dual_offset = 1.;
const double penalty_offset      = 10.;

static SLEQP_RETCODE
estimate_dual_values(SleqpTrialPointSolver* solver, SleqpIterate* iterate)
{
  SLEQP_CALL(sleqp_estimate_duals(solver->estimation_data,
                                  iterate,
                                  sleqp_iterate_cons_dual(iterate),
                                  sleqp_iterate_vars_dual(iterate)));

  SLEQP_CALL(sleqp_eqp_solver_add_violated_multipliers(solver->eqp_solver,
                                                       solver->multipliers));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
update_penalty(SleqpTrialPointSolver* solver)
{
  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;

  const int num_constraints = sleqp_problem_num_cons(problem);

  if (num_constraints == 0)
  {
    return SLEQP_OKAY;
  }

  const double feas_eps
    = sleqp_params_value(solver->params, SLEQP_PARAM_FEAS_TOL);

  bool is_feasible = sleqp_iterate_is_feasible(iterate,
                                               solver->feasibility_residuum,
                                               feas_eps);

  if (is_feasible)
  {
    solver->locally_infeasible = false;

    if (solver->allow_global_reset)
    {
      const SleqpVec* cons_dual = sleqp_iterate_cons_dual(iterate);
      const SleqpVec* vars_dual = sleqp_iterate_vars_dual(iterate);

      const double cons_dual_norm = sleqp_vec_inf_norm(cons_dual);
      const double vars_dual_norm = sleqp_vec_inf_norm(vars_dual);

      const double dual_norm = SLEQP_MAX(cons_dual_norm, vars_dual_norm);

      const double max_allowed_penalty
        = allowed_dual_factor * (dual_norm + allowed_dual_offset);

      if (solver->penalty_parameter > max_allowed_penalty)
      {
        sleqp_log_debug("Performing global penalty parameter reset");

        solver->penalty_parameter = dual_norm + penalty_offset;

        solver->performed_global_reset = true;
      }
    }
  }
  else
  {
    SLEQP_CALL(sleqp_update_penalty(problem,
                                    iterate,
                                    solver->cauchy_data,
                                    &(solver->penalty_parameter),
                                    &(solver->locally_infeasible)));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_cauchy_direction(SleqpTrialPointSolver* solver)
{
  SleqpIterate* iterate = solver->iterate;

  const bool enable_restoration
    = sleqp_options_bool_value(solver->options,
                               SLEQP_OPTION_BOOL_ENABLE_RESTORATION_PHASE);

  SLEQP_CALL(sleqp_cauchy_set_iterate(solver->cauchy_data,
                                      iterate,
                                      solver->lp_trust_radius));

  SLEQP_CALL(sleqp_cauchy_solve(solver->cauchy_data,
                                sleqp_iterate_obj_grad(iterate),
                                solver->penalty_parameter,
                                SLEQP_CAUCHY_OBJECTIVE_TYPE_DEFAULT));

  double criticality_bound;

  SLEQP_CALL(sleqp_cauchy_compute_criticality_bound(solver->cauchy_data,
                                                    solver->current_merit_value,
                                                    &criticality_bound));

  sleqp_log_debug("Criticality bound: %g", criticality_bound);

  SLEQP_CALL(
    sleqp_cauchy_get_direction(solver->cauchy_data, solver->cauchy_direction));

  SLEQP_CALL(sleqp_cauchy_get_working_set(solver->cauchy_data, iterate));

  const double original_penalty = solver->penalty_parameter;

  SLEQP_CALL(update_penalty(solver));

  if (solver->locally_infeasible && enable_restoration)
  {
    return SLEQP_OKAY;
  }

  if (original_penalty != solver->penalty_parameter)
  {
    SLEQP_CALL(sleqp_cauchy_get_direction(solver->cauchy_data,
                                          solver->cauchy_direction));

    SLEQP_CALL(sleqp_cauchy_get_working_set(solver->cauchy_data, iterate));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_cauchy_step_parametric(SleqpTrialPointSolver* solver,
                               double* cauchy_merit_value,
                               bool* full_step)
{
  SleqpIterate* iterate = solver->iterate;

  const double eps = sleqp_params_value(solver->params, SLEQP_PARAM_EPS);

  SLEQP_NUM_ASSERT_PARAM(eps);

  {
    SLEQP_CALL(compute_cauchy_direction(solver));

    SLEQP_CALL(sleqp_aug_jac_set_iterate(solver->aug_jac, iterate));

    SLEQP_CALL(sleqp_eqp_solver_set_iterate(solver->eqp_solver,
                                            iterate,
                                            solver->aug_jac,
                                            solver->trust_radius,
                                            solver->penalty_parameter));

    SLEQP_CALL(estimate_dual_values(solver, iterate));
  }

  {
    SLEQP_CALL(sleqp_parametric_solver_set_penalty(solver->parametric_solver,
                                                   solver->penalty_parameter));

    SLEQP_CALL(sleqp_parametric_solver_solve(solver->parametric_solver,
                                             iterate,
                                             solver->cauchy_data,
                                             solver->cauchy_step,
                                             solver->cauchy_hessian_step,
                                             solver->multipliers,
                                             &(solver->lp_trust_radius),
                                             cauchy_merit_value));
  }

  SLEQP_CALL(sleqp_working_set_copy(sleqp_iterate_working_set(iterate),
                                    solver->parametric_original_working_set));

  SLEQP_CALL(sleqp_cauchy_get_working_set(solver->cauchy_data, iterate));

  // Reconstruct the augmented Jacobian if required
  if (!sleqp_working_set_eq(solver->parametric_original_working_set,
                            sleqp_iterate_working_set(iterate)))
  {
    SLEQP_CALL(sleqp_aug_jac_set_iterate(solver->aug_jac, iterate));
  }

  SLEQP_CALL(sleqp_linesearch_set_iterate(solver->linesearch,
                                          iterate,
                                          solver->penalty_parameter,
                                          solver->trust_radius));

  SLEQP_CALL(sleqp_eqp_solver_set_iterate(solver->eqp_solver,
                                          iterate,
                                          solver->aug_jac,
                                          solver->trust_radius,
                                          solver->penalty_parameter));

#if !defined(NDEBUG)

  {
    double actual_quadratic_merit_value;

    double obj_dual = 1.;

    SLEQP_CALL(sleqp_merit_quadratic(solver->merit,
                                     iterate,
                                     &obj_dual,
                                     solver->cauchy_step,
                                     solver->multipliers,
                                     solver->penalty_parameter,
                                     &actual_quadratic_merit_value));

    sleqp_assert_is_eq(*cauchy_merit_value, actual_quadratic_merit_value, eps);
  }

#endif

  (*full_step) = true;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_cauchy_step_simple(SleqpTrialPointSolver* solver,
                           double* cauchy_merit_value,
                           bool quadratic_model,
                           bool* full_step)
{
  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;

  const double eps = sleqp_params_value(solver->params, SLEQP_PARAM_EPS);

  SLEQP_NUM_ASSERT_PARAM(eps);

  const double one = 1.;

  // compute Cauchy direction / step and dual estimation
  {
    SLEQP_CALL(compute_cauchy_direction(solver));

    SLEQP_CALL(sleqp_aug_jac_set_iterate(solver->aug_jac, iterate));

    SLEQP_CALL(sleqp_eqp_solver_set_iterate(solver->eqp_solver,
                                            iterate,
                                            solver->aug_jac,
                                            solver->trust_radius,
                                            solver->penalty_parameter));

    SLEQP_CALL(estimate_dual_values(solver, iterate));

#if !defined(NDEBUG)

    {
      bool in_working_set = false;

      SLEQP_CALL(sleqp_direction_in_working_set(problem,
                                                iterate,
                                                solver->cauchy_direction,
                                                solver->dense_cache,
                                                eps,
                                                &in_working_set));

      sleqp_num_assert(in_working_set);
    }

#endif

    SLEQP_CALL(sleqp_vec_copy(solver->cauchy_direction, solver->cauchy_step));

    if (!quadratic_model)
    {
      (*full_step) = true;

      return SLEQP_OKAY;
    }

    SLEQP_CALL(sleqp_problem_hess_prod(problem,
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
                                            full_step,
                                            cauchy_merit_value));

#if !defined(NDEBUG)

    {
      double actual_quadratic_merit_value, exact_iterate_value;

      double func_dual = 1.;

      SLEQP_CALL(sleqp_merit_quadratic(solver->merit,
                                       iterate,
                                       &func_dual,
                                       solver->cauchy_step,
                                       solver->multipliers,
                                       solver->penalty_parameter,
                                       &actual_quadratic_merit_value));

      sleqp_assert_is_eq(*cauchy_merit_value,
                         actual_quadratic_merit_value,
                         eps);

      SLEQP_CALL(sleqp_merit_func(solver->merit,
                                  iterate,
                                  solver->penalty_parameter,
                                  &exact_iterate_value));

      // quadratic merit at d = 0 corresponds to the
      // exact iterate value. The Cauchy step should
      // be at least as good
      sleqp_assert_is_leq(*cauchy_merit_value, exact_iterate_value, eps);
    }

#endif
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_trial_point_solver_compute_cauchy_step(SleqpTrialPointSolver* solver,
                                             double* cauchy_merit_value,
                                             bool quadratic_model,
                                             bool* full_step)
{
  SLEQP_PARAMETRIC_CAUCHY parametric_cauchy
    = sleqp_options_enum_value(solver->options,
                               SLEQP_OPTION_ENUM_PARAMETRIC_CAUCHY);

  SleqpTimer* timer = solver->elapsed_timer;
  double time_limit = solver->time_limit;

  double remaining_time = sleqp_timer_remaining_time(timer, time_limit);

  SLEQP_CALL(sleqp_cauchy_set_time_limit(solver->cauchy_data, remaining_time));

  if (parametric_cauchy != SLEQP_PARAMETRIC_CAUCHY_DISABLED)
  {
    SLEQP_CALL(
      compute_cauchy_step_parametric(solver, cauchy_merit_value, full_step));
  }
  else
  {
    SLEQP_CALL(compute_cauchy_step_simple(solver,
                                          cauchy_merit_value,
                                          quadratic_model,
                                          full_step));
  }

  return SLEQP_OKAY;
}
