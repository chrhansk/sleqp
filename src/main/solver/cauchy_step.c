#include "solver.h"

#include "cmp.h"
#include "fail.h"
#include "penalty.h"

static SLEQP_RETCODE estimate_dual_values(SleqpSolver* solver,
                                          SleqpIterate* iterate)
{
  SleqpOptions* options = solver->options;

  SLEQP_DUAL_ESTIMATION_TYPE estimation_type = sleqp_options_get_int(options,
                                                                     SLEQP_OPTION_INT_DUAL_ESTIMATION_TYPE);

  if(estimation_type == SLEQP_DUAL_ESTIMATION_TYPE_LSQ)
  {
    SLEQP_CALL(sleqp_dual_estimation_compute(solver->estimation_data,
                                             iterate,
                                             solver->estimation_residuals,
                                             solver->aug_jac));

#ifndef NDEBUG

    double unclipped_residuum = sleqp_sparse_vector_inf_norm(solver->estimation_residuals);
    double residuum;

    SLEQP_CALL(sleqp_iterate_stationarity_residuum(solver->problem,
                                                   solver->iterate,
                                                   solver->dense_cache,
                                                   &residuum));

    const double eps = sleqp_params_get(solver->params,
                                        SLEQP_PARAM_EPS);

    SLEQP_NUM_ASSERT_PARAM(eps);

    sleqp_assert_is_geq(residuum,
                        unclipped_residuum,
                        eps);

#endif

  }
  else
  {
    assert(estimation_type == SLEQP_DUAL_ESTIMATION_TYPE_LP);

    SLEQP_CALL(sleqp_cauchy_get_dual_estimation(solver->cauchy_data,
                                                iterate));
  }

  SLEQP_CALL(sleqp_eqp_solver_add_violated_multipliers(solver->eqp_solver,
                                                       solver->multipliers));

  return SLEQP_OKAY;
}


// Bound on the criticality measure used in
// "On the Convergence of Successive Linear Programming Algorithms"
static double compute_criticality_bound(SleqpSolver* solver)
{
  double objective_value;

  SLEQP_CALL(sleqp_cauchy_get_objective_value(solver->cauchy_data,
                                              &objective_value));

  const double reduction = solver->current_merit_value - objective_value;

  const double eps = sleqp_params_get(solver->params,
                                      SLEQP_PARAM_EPS);

  SLEQP_NUM_ASSERT_PARAM(eps);

  sleqp_assert_is_geq(reduction, 0., eps);

  const double criticality_bound = reduction / SLEQP_MIN(solver->lp_trust_radius, 1.);

  return criticality_bound;
}

static SLEQP_RETCODE
compute_cauchy_direction(SleqpSolver* solver)
{
  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;

  SLEQP_CALL(sleqp_cauchy_set_iterate(solver->cauchy_data,
                                      iterate,
                                      solver->lp_trust_radius));

  SLEQP_CALL(sleqp_cauchy_solve(solver->cauchy_data,
                                sleqp_iterate_get_func_grad(iterate),
                                solver->penalty_parameter,
                                SLEQP_CAUCHY_OBJECTIVE_TYPE_DEFAULT));

  const double criticality_bound = compute_criticality_bound(solver);

  sleqp_log_debug("Criticality bound: %g", criticality_bound);

  SLEQP_CALL(sleqp_cauchy_get_direction(solver->cauchy_data,
                                        solver->cauchy_direction));

  SLEQP_CALL(sleqp_cauchy_get_working_set(solver->cauchy_data,
                                          iterate));

  const double original_penalty = solver->penalty_parameter;

  SLEQP_CALL(sleqp_update_penalty(problem,
                                  iterate,
                                  solver->cauchy_data,
                                  &(solver->penalty_parameter),
                                  &(solver->locally_infeasible)));

  /*
  if(solver->locally_infeasible)
  {
    return SLEQP_OKAY;
  }
  */

  if(original_penalty != solver->penalty_parameter)
  {
    SLEQP_CALL(sleqp_cauchy_get_direction(solver->cauchy_data,
                                          solver->cauchy_direction));

    SLEQP_CALL(sleqp_cauchy_get_working_set(solver->cauchy_data,
                                            iterate));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_cauchy_step_parametric(SleqpSolver* solver,
                               double* cauchy_merit_value,
                               bool* full_step)
{
  SleqpIterate* iterate = solver->iterate;

  const double eps = sleqp_params_get(solver->params, SLEQP_PARAM_EPS);

  SLEQP_NUM_ASSERT_PARAM(eps);

  {
    SLEQP_CALL(compute_cauchy_direction(solver));

    SLEQP_CALL(sleqp_aug_jac_set_iterate(solver->aug_jac,
                                         iterate));

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

  SLEQP_CALL(sleqp_working_set_copy(sleqp_iterate_get_working_set(iterate),
                                    solver->parametric_original_working_set));

  SLEQP_CALL(sleqp_cauchy_get_working_set(solver->cauchy_data,
                                          iterate));

  // Reconstruct the augmented Jacobian if required
  if(!sleqp_working_set_eq(solver->parametric_original_working_set,
                           sleqp_iterate_get_working_set(iterate)))
  {
    SLEQP_CALL(sleqp_aug_jac_set_iterate(solver->aug_jac,
                                         iterate));
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
  }

#endif

  (*full_step) = true;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_cauchy_step_simple(SleqpSolver* solver,
                           double* cauchy_merit_value,
                           bool quadratic_model,
                           bool* full_step)
{
  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;

  const double eps = sleqp_params_get(solver->params,
                                      SLEQP_PARAM_EPS);

  SLEQP_NUM_ASSERT_PARAM(eps);

  const double zero_eps = sleqp_params_get(solver->params,
                                           SLEQP_PARAM_ZERO_EPS);

  const double one = 1.;

  // compute Cauchy direction / step and dual estimation
  {
    SLEQP_CALL(compute_cauchy_direction(solver));

    SLEQP_CALL(sleqp_aug_jac_set_iterate(solver->aug_jac,
                                         iterate));

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

    SLEQP_CALL(sleqp_sparse_vector_copy(solver->cauchy_direction,
                                        solver->cauchy_step));

    if(!quadratic_model)
    {
      (*full_step) = true;

      solver->cauchy_step_length = 1.;

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
                                            &solver->cauchy_step_length,
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
      sleqp_assert_is_leq(*cauchy_merit_value,
                          exact_iterate_value,
                          eps);
    }

#endif

    (*full_step) = sleqp_is_eq(solver->cauchy_step_length,
                               1.,
                               zero_eps);
  }

  return SLEQP_OKAY;
}


SLEQP_RETCODE sleqp_solver_compute_cauchy_step(SleqpSolver* solver,
                                               double* cauchy_merit_value,
                                               bool quadratic_model,
                                               bool* full_step)
{
  SLEQP_PARAMETRIC_CAUCHY parametric_cauchy = sleqp_options_get_int(solver->options,
                                                                    SLEQP_OPTION_INT_PARAMETRIC_CAUCHY);

  if(parametric_cauchy != SLEQP_PARAMETRIC_CAUCHY_DISABLED)
  {
    SLEQP_CALL(compute_cauchy_step_parametric(solver,
                                              cauchy_merit_value,
                                              full_step));
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
