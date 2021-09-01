#include "solver.h"

#include "dyn.h"
#include "fail.h"

static SLEQP_RETCODE
compute_trial_iterate_from_direction(SleqpSolver* solver,
                                     SleqpSparseVec* direction)
{
  SleqpProblem* problem = solver->problem;

  SleqpIterate* iterate = solver->iterate;
  SleqpIterate* trial_iterate = solver->trial_iterate;

  const double zero_eps = sleqp_params_get(solver->params,
                                           SLEQP_PARAM_ZERO_EPS);

  SLEQP_CALL(sleqp_sparse_vector_add(sleqp_iterate_get_primal(iterate),
                                     direction,
                                     zero_eps,
                                     solver->initial_trial_point));

  SLEQP_CALL(sleqp_sparse_vector_clip(solver->initial_trial_point,
                                      sleqp_problem_var_lb(problem),
                                      sleqp_problem_var_ub(problem),
                                      zero_eps,
                                      sleqp_iterate_get_primal(trial_iterate)));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_trial_point_deterministic(SleqpSolver* solver,
                                  double* trial_merit_value,
                                  bool* full_step)
{
  const SleqpOptions* options = solver->options;

  const bool quadratic_model = sleqp_options_get_bool(options,
                                                      SLEQP_OPTION_BOOL_USE_QUADRATIC_MODEL);

  const bool perform_newton_step = quadratic_model &&
    sleqp_options_get_bool(options, SLEQP_OPTION_BOOL_PERFORM_NEWTON_STEP);

  if(perform_newton_step)
  {
    SLEQP_CALL(sleqp_solver_compute_trial_point_newton(solver,
                                                       trial_merit_value,
                                                       full_step));
  }
  else
  {
    SLEQP_CALL(sleqp_solver_compute_trial_point_simple(solver,
                                                       trial_merit_value,
                                                       quadratic_model,
                                                       full_step));
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solver_compute_trial_point_simple(SleqpSolver* solver,
                                                      double* cauchy_merit_value,
                                                      bool quadratic_model,
                                                      bool* full_step)
{
  SleqpIterate* iterate = solver->iterate;

  const double eps = sleqp_params_get(solver->params,
                                      SLEQP_PARAM_EPS);

  SLEQP_NUM_ASSERT_PARAM(eps);

  SLEQP_CALL(sleqp_solver_compute_cauchy_step(solver,
                                              cauchy_merit_value,
                                              quadratic_model,
                                              full_step));

  const SleqpSparseVec* trial_step = solver->cauchy_step;

  // Compute merit value
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

  SLEQP_CALL(compute_trial_iterate_from_direction(solver,
                                                  solver->cauchy_step));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solver_compute_trial_point_newton(SleqpSolver* solver,
                                                      double* trial_merit_value,
                                                      bool* full_step)
{
  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;

  const double eps = sleqp_params_get(solver->params, SLEQP_PARAM_EPS);

  SLEQP_NUM_ASSERT_PARAM(eps);

  const double one = 1.;

  double cauchy_merit_value;

  SLEQP_CALL(sleqp_solver_compute_cauchy_step(solver,
                                              &cauchy_merit_value,
                                              true,
                                              full_step));

  SLEQP_TR_SOLVER tr_solver = sleqp_options_get_int(solver->options,
                                                    SLEQP_OPTION_INT_TR_SOLVER);

  if(tr_solver == SLEQP_TR_SOLVER_LSQR)
  {
    SLEQP_CALL(sleqp_lsqr_solver_compute_step(solver->lsqr_solver,
                                              solver->newton_step));
  }
  else
  // compute Newton step
  {
    SLEQP_CALL(sleqp_newton_set_time_limit(solver->newton_data,
                                           sleqp_solver_remaining_time(solver)));

    SLEQP_CALL(sleqp_newton_compute_step(solver->newton_data,
                                         solver->multipliers,
                                         solver->newton_step));
  }

  SLEQP_CALL(sleqp_problem_hess_prod(problem,
                                     &one,
                                     solver->newton_step,
                                     solver->multipliers,
                                     solver->newton_hessian_step));

  {
    SLEQP_LINESEARCH lineserach = sleqp_options_get_int(solver->options,
                                                        SLEQP_OPTION_INT_LINESEARCH);

    double step_length;

    if(lineserach == SLEQP_LINESEARCH_EXACT)
    {
      SLEQP_CALL(sleqp_linesearch_trial_step_exact(solver->linesearch,
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
    else
    {
      assert(lineserach == SLEQP_LINESEARCH_APPROX);

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

  SLEQP_CALL(compute_trial_iterate_from_direction(solver,
                                                  solver->trial_step));

  return SLEQP_OKAY;
}

static double compute_required_accuracy(SleqpSolver* solver,
                                        double model_reduction)
{
  const double accepted_reduction = sleqp_params_get(solver->params,
                                                     SLEQP_PARAM_ACCEPTED_REDUCTION);

  // TODO: Make this adjustable
  // must be > 0, < .5 *accepted_reduction
  const double required_accuracy_factor = .4 * accepted_reduction;

  return required_accuracy_factor * model_reduction;
}

static
SLEQP_RETCODE solver_refine_step(SleqpSolver* solver,
                                 double* model_trial_value,
                                 bool* full_step)
{
  SleqpProblem* problem = solver->problem;
  SleqpFunc* func = sleqp_problem_func(problem);

  SleqpOptions* options = solver->options;

  SleqpIterate* iterate = solver->iterate;

  assert(sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_DYNAMIC);

  const bool quadratic_model = sleqp_options_get_bool(options,
                                                      SLEQP_OPTION_BOOL_USE_QUADRATIC_MODEL);

  const bool perform_newton_step = quadratic_model &&
    sleqp_options_get_bool(options, SLEQP_OPTION_BOOL_PERFORM_NEWTON_STEP);

  while(true)
  {
    double current_accuracy;

    SLEQP_CALL(sleqp_dyn_func_get_accuracy(func, &current_accuracy));

    const double model_reduction = solver->current_merit_value - (*model_trial_value);

    const double required_accuracy = compute_required_accuracy(solver, model_reduction);

    if(current_accuracy <= required_accuracy)
    {
      break;
    }

    sleqp_log_debug("Current accuracy of %e is insufficient, reducing to %e",
                    current_accuracy,
                    required_accuracy);

    SLEQP_CALL(sleqp_dyn_func_set_accuracy(func, required_accuracy));

    // TODO: We do not really need to *set* the function value again...
    SLEQP_CALL(sleqp_set_and_evaluate(problem,
                                      iterate,
                                      SLEQP_VALUE_REASON_INIT));

    SLEQP_CALL(sleqp_merit_func(solver->merit_data,
                                iterate,
                                solver->penalty_parameter,
                                &solver->current_merit_value));

    // TODO: recompute more efficiently if possible
    if(perform_newton_step)
    {
      SLEQP_CALL(sleqp_solver_compute_trial_point_newton(solver,
                                                         model_trial_value,
                                                         full_step));
    }
    else
    {
      SLEQP_CALL(sleqp_solver_compute_trial_point_simple(solver,
                                                         model_trial_value,
                                                         quadratic_model,
                                                         full_step));
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solver_compute_trial_point_dynamic(SleqpSolver* solver,
                                                       double* trial_merit_value,
                                                       bool* full_step)
{
  SLEQP_CALL(compute_trial_point_deterministic(solver,
                                               trial_merit_value,
                                               full_step));

  SLEQP_CALL(solver_refine_step(solver,
                                trial_merit_value,
                                full_step));

  SLEQP_CALL(compute_trial_iterate_from_direction(solver,
                                                  solver->trial_step));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solver_compute_trial_point(SleqpSolver* solver,
                                               double* trial_merit_value,
                                               bool* full_step,
                                               bool* reject)
{
  SleqpProblem* problem = solver->problem;

  SleqpFunc* func = sleqp_problem_func(problem);

  // TODO: enable manual rejects
  *reject = false;

  if(sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_DYNAMIC)
  {
    SLEQP_CALL(sleqp_solver_compute_trial_point_dynamic(solver,
                                                        trial_merit_value,
                                                        full_step));

    return SLEQP_OKAY;
  }
  else
  {
    SLEQP_CALL(compute_trial_point_deterministic(solver,
                                                 trial_merit_value,
                                                 full_step));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_trial_point_soc_deterministic(SleqpSolver* solver,
                                      bool* reject)
{
  SleqpIterate* iterate = solver->iterate;
  SleqpIterate* trial_iterate = solver->trial_iterate;

  *reject = false;

  SLEQP_CALL(sleqp_soc_compute_step(solver->soc_data,
                                    solver->aug_jacobian,
                                    iterate,
                                    solver->trial_step,
                                    trial_iterate,
                                    solver->soc_step));

  SLEQP_CALL(compute_trial_iterate_from_direction(solver,
                                                  solver->soc_step));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_trial_point_soc_dynamic(SleqpSolver* solver,
                                bool* reject)
{
  SleqpProblem* problem = solver->problem;
  SleqpFunc* func = sleqp_problem_func(problem);

  assert(sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_DYNAMIC);

  SleqpIterate* iterate = solver->iterate;
  SleqpIterate* trial_iterate = solver->trial_iterate;

  SleqpOptions* options = solver->options;

  const bool quadratic_model = sleqp_options_get_bool(options,
                                                      SLEQP_OPTION_BOOL_USE_QUADRATIC_MODEL);

  SLEQP_CALL(sleqp_soc_compute_step(solver->soc_data,
                                    solver->aug_jacobian,
                                    iterate,
                                    solver->trial_step,
                                    trial_iterate,
                                    solver->soc_step));

  double soc_model_merit = SLEQP_NONE;

  if(quadratic_model)
  {
    const double one = 1.;

    SLEQP_CALL(sleqp_merit_quadratic(solver->merit_data,
                                     iterate,
                                     &one,
                                     solver->soc_step,
                                     solver->multipliers,
                                     solver->penalty_parameter,
                                     &soc_model_merit));
  }
  else
  {
    SLEQP_CALL(sleqp_merit_linear(solver->merit_data,
                                  iterate,
                                  solver->soc_step,
                                  solver->penalty_parameter,
                                  &soc_model_merit));
  }

  const double model_reduction = solver->current_merit_value - soc_model_merit;

  double current_accuracy;

  SLEQP_CALL(sleqp_dyn_func_get_accuracy(func, &current_accuracy));

  const double required_accuracy = compute_required_accuracy(solver,
                                                             model_reduction);

  *reject = (current_accuracy > required_accuracy);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solver_compute_trial_point_soc(SleqpSolver* solver,
                                                   bool* reject)
{
  SleqpProblem* problem = solver->problem;
  SleqpFunc* func = sleqp_problem_func(problem);

  if(sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_DYNAMIC)
  {
    SLEQP_CALL(compute_trial_point_soc_dynamic(solver,
                                               reject));
  }
  else
  {
    SLEQP_CALL(compute_trial_point_soc_deterministic(solver,
                                                     reject));
  }

  return SLEQP_OKAY;
}
