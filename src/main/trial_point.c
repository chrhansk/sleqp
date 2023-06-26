#include "trial_point.h"

#include "direction.h"
#include "dyn.h"
#include "error.h"
#include "fail.h"
#include "func.h"
#include "gauss_newton.h"
#include "mem.h"
#include "newton.h"
#include "settings.h"
#include "soc.h"

#include "aug_jac/box_constrained_aug_jac.h"
#include "aug_jac/direct_aug_jac.h"
#include "aug_jac/reduced_aug_jac.h"
#include "aug_jac/standard_aug_jac.h"
#include "aug_jac/unconstrained_aug_jac.h"

#include "cauchy/box_constrained_cauchy.h"
#include "cauchy/standard_cauchy.h"
#include "cauchy/unconstrained_cauchy.h"

#include "dual_estimation/dual_estimation_lp.h"
#include "dual_estimation/dual_estimation_lsq.h"
#include "dual_estimation/dual_estimation_mixed.h"

#include "fact/fact.h"
#include "fact/fact_qr.h"

static SLEQP_RETCODE
create_dual_estimation(SleqpTrialPointSolver* solver)
{
  SleqpSettings* settings = solver->settings;

  SLEQP_DUAL_ESTIMATION_TYPE estimation_type
    = sleqp_settings_enum_value(settings,
                                SLEQP_SETTINGS_ENUM_DUAL_ESTIMATION_TYPE);

  if (estimation_type == SLEQP_DUAL_ESTIMATION_TYPE_LP)
  {
    SLEQP_CALL(sleqp_dual_estimation_lp_create(&solver->estimation_data,
                                               solver->cauchy_data));
  }
  else if (estimation_type == SLEQP_DUAL_ESTIMATION_TYPE_LSQ)
  {
    SLEQP_CALL(sleqp_dual_estimation_lsq_create(&solver->estimation_data,
                                                solver->problem,
                                                solver->aug_jac));
  }
  else
  {
    assert(estimation_type == SLEQP_DUAL_ESTIMATION_TYPE_MIXED);

    SLEQP_CALL(sleqp_dual_estimation_mixed_create(&solver->estimation_data,
                                                  solver->problem,
                                                  solver->cauchy_data,
                                                  solver->aug_jac));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
create_aug_jac(SleqpTrialPointSolver* solver)
{
  SleqpProblem* problem   = solver->problem;
  SleqpSettings* settings = solver->settings;

  const int num_constraints = sleqp_problem_num_cons(problem);

  if (sleqp_problem_is_unconstrained(problem))
  {
    SLEQP_CALL(sleqp_unconstrained_aug_jac_create(&solver->aug_jac, problem));
  }
  else if (num_constraints == 0)
  {
    SLEQP_CALL(sleqp_box_constrained_aug_jac_create(&solver->aug_jac, problem));
  }
  else
  {
    // create sparse factorization

    const SLEQP_AUG_JAC_METHOD aug_jac_method
      = sleqp_settings_enum_value(settings, SLEQP_SETTINGS_ENUM_AUG_JAC_METHOD);

    bool requires_psd = false;

    switch (aug_jac_method)
    {
    case SLEQP_AUG_JAC_AUTO:
      SLEQP_CALL(sleqp_fact_create_default(&solver->fact, settings));

      requires_psd = (sleqp_fact_flags(solver->fact) & SLEQP_FACT_FLAGS_PSD);

      if (requires_psd)
      {
        SLEQP_CALL(sleqp_reduced_aug_jac_create(&solver->aug_jac,
                                                problem,
                                                settings,
                                                solver->fact));
      }
      else
      {
        SLEQP_CALL(sleqp_standard_aug_jac_create(&solver->aug_jac,
                                                 problem,
                                                 settings,
                                                 solver->fact));
      }

      break;
    case SLEQP_AUG_JAC_STANDARD:
      SLEQP_CALL(sleqp_fact_create_default(&solver->fact, settings));

      requires_psd = (sleqp_fact_flags(solver->fact) & SLEQP_FACT_FLAGS_PSD);

      if (requires_psd)
      {
        sleqp_raise(SLEQP_ILLEGAL_ARGUMENT,
                    "Factorization requires PSD matrices");
      }

      SLEQP_CALL(sleqp_standard_aug_jac_create(&solver->aug_jac,
                                               problem,
                                               settings,
                                               solver->fact));
      break;
    case SLEQP_AUG_JAC_REDUCED:
      SLEQP_CALL(sleqp_fact_create_default(&solver->fact, settings));

      SLEQP_CALL(sleqp_reduced_aug_jac_create(&solver->aug_jac,
                                              problem,
                                              settings,
                                              solver->fact));

      break;
    case SLEQP_AUG_JAC_DIRECT:
#ifdef SLEQP_HAVE_QR_FACT

      ;
      SleqpFactQR* qr_fact = NULL;

      SLEQP_CALL(sleqp_fact_qr_create_default(&qr_fact, settings));

      SLEQP_CALL(sleqp_direct_aug_jac_create(&solver->aug_jac,
                                             problem,
                                             settings,
                                             qr_fact));

      SLEQP_CALL(sleqp_qr_release(&qr_fact));
#else
      sleqp_raise(SLEQP_ILLEGAL_ARGUMENT, "No QR factorization available");
#endif
      break;
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
create_cauchy_solver(SleqpTrialPointSolver* solver)
{
  SleqpProblem* problem   = solver->problem;
  SleqpSettings* settings = solver->settings;

  const int num_constraints = sleqp_problem_num_cons(problem);

  if (sleqp_problem_is_unconstrained(problem))
  {
    SLEQP_CALL(sleqp_unconstrained_cauchy_create(&solver->cauchy_data,
                                                 problem,
                                                 settings));
  }
  else if (num_constraints == 0)
  {
    SLEQP_CALL(sleqp_box_constrained_cauchy_create(&solver->cauchy_data,
                                                   problem,
                                                   settings));
  }
  else
  {
    SLEQP_CALL(
      sleqp_standard_cauchy_create(&solver->cauchy_data, problem, settings));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
create_eqp_solver(SleqpTrialPointSolver* solver)
{
  SleqpProblem* problem   = solver->problem;
  SleqpSettings* settings = solver->settings;

  SLEQP_TR_SOLVER tr_solver
    = sleqp_settings_enum_value(settings, SLEQP_SETTINGS_ENUM_TR_SOLVER);

  if (tr_solver == SLEQP_TR_SOLVER_LSQR)
  {
    SleqpFunc* func = sleqp_problem_func(problem);

    if (sleqp_func_get_type(func) != SLEQP_FUNC_TYPE_LSQ)
    {
      sleqp_raise(SLEQP_ILLEGAL_ARGUMENT,
                  "LSQR solver is only available for LSQ problems");
    }

    SLEQP_CALL(sleqp_gauss_newton_solver_create(&solver->eqp_solver,
                                                solver->problem,
                                                settings,
                                                solver->working_step));
  }
  else
  {
    SLEQP_CALL(sleqp_newton_solver_create(&solver->eqp_solver,
                                          solver->problem,
                                          settings,
                                          solver->working_step));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
create_parametric_solver(SleqpTrialPointSolver* solver)
{
  SLEQP_PARAMETRIC_CAUCHY parametric_cauchy
    = sleqp_settings_enum_value(solver->settings,
                                SLEQP_SETTINGS_ENUM_PARAMETRIC_CAUCHY);

  if (parametric_cauchy == SLEQP_PARAMETRIC_CAUCHY_DISABLED)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_parametric_solver_create(&solver->parametric_solver,
                                            solver->problem,
                                            solver->settings,
                                            solver->merit,
                                            solver->linesearch));

  SLEQP_CALL(sleqp_working_set_create(&solver->parametric_original_working_set,
                                      solver->problem));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_trial_point_solver_create(SleqpTrialPointSolver** star,
                                SleqpProblem* problem,
                                SleqpSettings* settings)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpTrialPointSolver* solver = *star;

  *solver = (SleqpTrialPointSolver){0};

  solver->refcount = 1;

  SLEQP_CALL(sleqp_problem_capture(problem));
  solver->problem = problem;

  SLEQP_CALL(sleqp_settings_capture(settings));
  solver->settings = settings;

  const int num_variables   = sleqp_problem_num_vars(solver->problem);
  const int num_constraints = sleqp_problem_num_cons(solver->problem);

  SLEQP_CALL(sleqp_vec_create_empty(&solver->lp_step, num_variables));

  SLEQP_CALL(
    sleqp_direction_create(&solver->cauchy_direction, problem, settings));

  SLEQP_CALL(
    sleqp_vec_create_empty(&solver->estimation_residuals, num_variables));

  SLEQP_CALL(
    sleqp_direction_create(&solver->newton_direction, problem, settings));

  SLEQP_CALL(sleqp_direction_create(&solver->soc_direction, problem, settings));

  SLEQP_CALL(
    sleqp_direction_create(&solver->trial_direction, problem, settings));

  SLEQP_CALL(sleqp_vec_create_empty(&solver->multipliers, num_constraints));

  SLEQP_CALL(
    sleqp_vec_create_empty(&solver->initial_trial_point, num_variables));

  SLEQP_CALL(sleqp_merit_create(&solver->merit, problem, settings));

  SLEQP_CALL(create_cauchy_solver(solver));

  SLEQP_CALL(create_aug_jac(solver));

  SLEQP_CALL(create_dual_estimation(solver));

  SLEQP_CALL(sleqp_linesearch_create(&solver->linesearch,
                                     solver->problem,
                                     settings,
                                     solver->merit));

  SLEQP_CALL(sleqp_working_step_create(&solver->working_step,
                                       solver->problem,
                                       settings));

  SLEQP_CALL(create_eqp_solver(solver));

  SLEQP_CALL(
    sleqp_soc_data_create(&solver->soc_data, solver->problem, settings));

  SLEQP_CALL(create_parametric_solver(solver));

  SLEQP_CALL(sleqp_alloc_array(&solver->dense_cache,
                               SLEQP_MAX(num_variables, num_constraints)));

  SLEQP_CALL(sleqp_timer_create(&solver->elapsed_timer));

  solver->time_limit = SLEQP_NONE;

  solver->penalty_parameter = SLEQP_NONE;
  solver->trust_radius      = SLEQP_NONE;
  solver->lp_trust_radius   = SLEQP_NONE;

  {
    SleqpFunc* func = sleqp_problem_func(problem);

    if (sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_DYNAMIC)
    {
      SLEQP_CALL(sleqp_alloc_array(&solver->cons_weights, num_constraints));
      SLEQP_CALL(sleqp_dyn_func_set_error_bound(func, 1.));
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_trial_point_solver_set_iterate(SleqpTrialPointSolver* solver,
                                     SleqpIterate* iterate)
{
  SLEQP_CALL(sleqp_iterate_release(&solver->iterate));

  SLEQP_CALL(sleqp_iterate_capture(iterate));
  solver->iterate = iterate;

  SLEQP_CALL(sleqp_merit_func(solver->merit,
                              iterate,
                              solver->penalty_parameter,
                              &solver->current_merit_value));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_trial_point_solver_set_penalty_info(SleqpTrialPointSolver* solver,
                                          double feas_res,
                                          bool allow_global_reset)
{
  solver->feasibility_residuum   = feas_res;
  solver->allow_global_reset     = allow_global_reset;
  solver->performed_global_reset = false;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_trial_point_solver_set_trust_radius(SleqpTrialPointSolver* solver,
                                          double trust_radius)
{
  assert(trust_radius > 0.);

  solver->trust_radius = trust_radius;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_trial_point_solver_set_lp_trust_radius(SleqpTrialPointSolver* solver,
                                             double lp_trust_radius)
{
  assert(lp_trust_radius > 0.);

  solver->lp_trust_radius = lp_trust_radius;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_trial_point_solver_set_penalty(SleqpTrialPointSolver* solver,
                                     double penalty_parameter)
{
  assert(penalty_parameter > 0.);

  const double last_penalty_parameter = solver->penalty_parameter;

  solver->penalty_parameter = penalty_parameter;

  if (solver->penalty_parameter != last_penalty_parameter)
  {
    SLEQP_CALL(sleqp_trial_point_solver_set_cons_weights(solver));
  }

  SleqpIterate* iterate = solver->iterate;

  SLEQP_CALL(sleqp_merit_func(solver->merit,
                              iterate,
                              solver->penalty_parameter,
                              &solver->current_merit_value));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_trial_point_solver_penalty(SleqpTrialPointSolver* solver,
                                 double* penalty_parameter)
{
  *penalty_parameter = solver->penalty_parameter;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_trial_point_solver_set_cons_weights(SleqpTrialPointSolver* solver)
{
  SleqpProblem* problem = solver->problem;
  SleqpFunc* func       = sleqp_problem_func(problem);

  assert(solver->penalty_parameter != SLEQP_NONE);

  SLEQP_CALL(sleqp_dyn_set_penalty_cons_weights(func,
                                                solver->penalty_parameter,
                                                solver->cons_weights));

  return SLEQP_OKAY;
}

bool
sleqp_trial_point_solver_locally_infeasible(SleqpTrialPointSolver* solver)
{
  SleqpSettings* settings = solver->settings;

  const double feas_tol
    = sleqp_settings_real_value(settings, SLEQP_SETTINGS_REAL_FEAS_TOL);
  const double eps
    = sleqp_settings_real_value(settings, SLEQP_SETTINGS_REAL_EPS);

  // Iterate is feasible anyways
  if (solver->feasibility_residuum <= feas_tol)
  {
    return false;
  }

  const double lp_step_norm = sleqp_vec_norm(solver->lp_step);

  // Nonzero step computed
  if (!sleqp_is_zero(lp_step_norm, eps))
  {
    return false;
  }

  const double step_norm
    = sleqp_vec_norm(sleqp_direction_primal(solver->trial_direction));

  // Step may be zero, but may be possible to
  // escape using second order information
  if (!sleqp_is_zero(step_norm, eps))
  {
    return false;
  }

  return true;
}

SLEQP_RETCODE
sleqp_trial_point_solver_penalty_info(SleqpTrialPointSolver* solver,
                                      bool* performed_global_reset)
{
  *performed_global_reset = solver->performed_global_reset;
  return SLEQP_OKAY;
}

SleqpVec*
sleqp_trial_point_solver_multipliers(SleqpTrialPointSolver* solver)
{
  return solver->multipliers;
}

SleqpVec*
sleqp_trial_point_solver_cauchy_step(SleqpTrialPointSolver* solver)
{
  return sleqp_direction_primal(solver->cauchy_direction);
}

SleqpDirection*
sleqp_trial_point_solver_trial_direction(SleqpTrialPointSolver* solver)
{
  return solver->trial_direction;
}

SleqpVec*
sleqp_trial_point_solver_soc_step(SleqpTrialPointSolver* solver)
{
  return sleqp_direction_primal(solver->soc_direction);
}

SLEQP_RETCODE
sleqp_trial_point_solver_rayleigh(SleqpTrialPointSolver* solver,
                                  double* min_rayleigh,
                                  double* max_rayleigh)
{
  SLEQP_CALL(sleqp_eqp_solver_current_rayleigh(solver->eqp_solver,
                                               min_rayleigh,
                                               max_rayleigh));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_trial_iterate_from_step(SleqpTrialPointSolver* solver,
                                const SleqpVec* step,
                                SleqpIterate* trial_iterate)
{
  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;

  const double zero_eps
    = sleqp_settings_real_value(solver->settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  SLEQP_CALL(sleqp_vec_add(sleqp_iterate_primal(iterate),
                           step,
                           zero_eps,
                           solver->initial_trial_point));

  SLEQP_CALL(sleqp_vec_clip(solver->initial_trial_point,
                            sleqp_problem_vars_lb(problem),
                            sleqp_problem_vars_ub(problem),
                            zero_eps,
                            sleqp_iterate_primal(trial_iterate)));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_trial_point_simple(SleqpTrialPointSolver* solver,
                           SleqpIterate* trial_iterate,
                           double* cauchy_merit_value,
                           bool quadratic_model,
                           bool* full_step)
{
  SLEQP_CALL(sleqp_trial_point_solver_compute_cauchy_step(solver,
                                                          cauchy_merit_value,
                                                          quadratic_model,
                                                          full_step));

  SleqpDirection* trial_direction = solver->cauchy_direction;

  // Compute merit value
  {
    SLEQP_CALL(sleqp_merit_linear(solver->merit,
                                  solver->iterate,
                                  trial_direction,
                                  solver->penalty_parameter,
                                  cauchy_merit_value));
  }

  if (quadratic_model)
  {
    double hessian_prod;

    SLEQP_CALL(sleqp_vec_dot(sleqp_direction_primal(trial_direction),
                             sleqp_direction_hess(trial_direction),
                             &hessian_prod));

    (*cauchy_merit_value) += .5 * hessian_prod;
  }

  SLEQP_CALL(
    sleqp_direction_copy(solver->cauchy_direction, solver->trial_direction));

  SLEQP_CALL(compute_trial_iterate_from_step(
    solver,
    sleqp_direction_primal(solver->trial_direction),
    trial_iterate));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_trial_point_newton(SleqpTrialPointSolver* solver,
                           SleqpIterate* trial_iterate,
                           double* trial_merit_value,
                           bool* failed_eqp_step,
                           bool* full_step)
{
  SleqpIterate* iterate = solver->iterate;

  SleqpTimer* timer = solver->elapsed_timer;
  double time_limit = solver->time_limit;

  double remaining_time = sleqp_timer_remaining_time(timer, time_limit);

  double cauchy_merit_value;

  SLEQP_CALL(sleqp_trial_point_solver_compute_cauchy_step(solver,
                                                          &cauchy_merit_value,
                                                          true,
                                                          full_step));

  SLEQP_CALL(
    sleqp_eqp_solver_set_time_limit(solver->eqp_solver, remaining_time));

  SLEQP_CALL(sleqp_eqp_solver_compute_direction(solver->eqp_solver,
                                                solver->multipliers,
                                                solver->newton_direction));

#if SLEQP_DEBUG

  {
    bool direction_valid;

    const double zero_eps
      = sleqp_settings_real_value(solver->settings,
                                  SLEQP_SETTINGS_REAL_ZERO_EPS);

    SLEQP_CALL(sleqp_direction_check(solver->cauchy_direction,
                                     solver->problem,
                                     iterate,
                                     solver->multipliers,
                                     solver->dense_cache,
                                     zero_eps,
                                     &direction_valid));

    sleqp_num_assert(direction_valid);
  }

  {
    bool direction_valid;

    const double zero_eps
      = sleqp_settings_real_value(solver->settings,
                                  SLEQP_SETTINGS_REAL_ZERO_EPS);

    SLEQP_CALL(sleqp_direction_check(solver->newton_direction,
                                     solver->problem,
                                     iterate,
                                     solver->multipliers,
                                     solver->dense_cache,
                                     zero_eps,
                                     &direction_valid));

    sleqp_num_assert(direction_valid);
  }

#endif

  {
    SLEQP_LINESEARCH lineserach
      = sleqp_settings_enum_value(solver->settings,
                                  SLEQP_SETTINGS_ENUM_LINESEARCH);

    double step_length;

    if (lineserach == SLEQP_LINESEARCH_EXACT)
    {
      SLEQP_CALL(sleqp_linesearch_trial_step_exact(solver->linesearch,
                                                   solver->cauchy_direction,
                                                   cauchy_merit_value,
                                                   solver->newton_direction,
                                                   solver->multipliers,
                                                   solver->trial_direction,
                                                   &step_length,
                                                   trial_merit_value));
    }
    else
    {
      assert(lineserach == SLEQP_LINESEARCH_APPROX);

      SLEQP_CALL(sleqp_linesearch_trial_step(solver->linesearch,
                                             solver->cauchy_direction,
                                             cauchy_merit_value,
                                             solver->newton_direction,
                                             solver->multipliers,
                                             solver->trial_direction,
                                             &step_length,
                                             trial_merit_value));
    }

    *(failed_eqp_step) = (step_length == 0.);
  }

#if SLEQP_DEBUG
  if (*failed_eqp_step)
  {
    assert(*trial_merit_value == cauchy_merit_value);
  }
#endif

  SLEQP_CALL(compute_trial_iterate_from_step(
    solver,
    sleqp_direction_primal(solver->trial_direction),
    trial_iterate));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_trial_point_deterministic(SleqpTrialPointSolver* solver,
                                  SleqpIterate* trial_iterate,
                                  double* trial_merit_value,
                                  bool* failed_eqp_step,
                                  bool* full_step)
{
  const SleqpSettings* settings = solver->settings;

  const bool quadratic_model
    = sleqp_settings_bool_value(settings,
                                SLEQP_SETTINGS_BOOL_USE_QUADRATIC_MODEL);

  const bool perform_newton_step
    = quadratic_model
      && sleqp_settings_bool_value(settings,
                                   SLEQP_SETTINGS_BOOL_PERFORM_NEWTON_STEP);

  if (perform_newton_step)
  {
    SLEQP_CALL(compute_trial_point_newton(solver,
                                          trial_iterate,
                                          trial_merit_value,
                                          failed_eqp_step,
                                          full_step));
  }
  else
  {
    SLEQP_CALL(compute_trial_point_simple(solver,
                                          trial_iterate,
                                          trial_merit_value,
                                          quadratic_model,
                                          full_step));
  }

#if SLEQP_DEBUG

  {
    bool direction_valid;

    const double zero_eps
      = sleqp_settings_real_value(solver->settings,
                                  SLEQP_SETTINGS_REAL_ZERO_EPS);

    SLEQP_CALL(sleqp_direction_check(solver->trial_direction,
                                     solver->problem,
                                     solver->iterate,
                                     solver->multipliers,
                                     solver->dense_cache,
                                     zero_eps,
                                     &direction_valid));

    sleqp_num_assert(direction_valid);
  }

  {
    double actual_merit_value = 0.;

    if (quadratic_model)
    {
      SLEQP_CALL(sleqp_merit_quadratic(solver->merit,
                                       solver->iterate,
                                       solver->trial_direction,
                                       solver->penalty_parameter,
                                       &actual_merit_value));
    }
    else
    {
      SLEQP_CALL(sleqp_merit_linear(solver->merit,
                                    solver->iterate,
                                    solver->trial_direction,
                                    solver->penalty_parameter,
                                    &actual_merit_value));
    }

    const double eps
      = sleqp_settings_real_value(solver->settings, SLEQP_SETTINGS_REAL_EPS);

    SLEQP_NUM_ASSERT_PARAM(eps);

    sleqp_assert_is_eq(*trial_merit_value, actual_merit_value, eps);
  }

#endif

  return SLEQP_OKAY;
}

static double
compute_required_error_bound(SleqpTrialPointSolver* solver,
                             double model_reduction)
{
  const double accepted_reduction
    = sleqp_settings_real_value(solver->settings,
                                SLEQP_SETTINGS_REAL_ACCEPTED_REDUCTION);

  // TODO: Make this adjustable
  // must be > 0, < .5 *accepted_reduction
  const double required_accuracy_factor = .4 * accepted_reduction;

  return required_accuracy_factor * model_reduction;
}

static SLEQP_RETCODE
refine_iterate(SleqpTrialPointSolver* solver,
               SleqpProblem* problem,
               SleqpIterate* iterate,
               double required_accuracy)
{
  SleqpFunc* func = sleqp_problem_func(problem);

  assert(sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_DYNAMIC);

  SLEQP_CALL(sleqp_dyn_func_set_error_bound(func, required_accuracy));

  double obj_val;

  SleqpVec* obj_grad = sleqp_iterate_obj_grad(iterate);
  SleqpMat* cons_jac = sleqp_iterate_cons_jac(iterate);
  SleqpVec* cons_val = sleqp_iterate_cons_val(iterate);

  SLEQP_CALL(
    sleqp_problem_eval(problem, &obj_val, obj_grad, cons_val, cons_jac));

  SLEQP_CALL(sleqp_iterate_set_obj_val(iterate, obj_val));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
solver_refine_step(SleqpTrialPointSolver* solver,
                   SleqpIterate* trial_iterate,
                   double* model_trial_value,
                   bool* full_step)
{
  SleqpProblem* problem = solver->problem;
  SleqpFunc* func       = sleqp_problem_func(problem);

  SleqpSettings* settings = solver->settings;

  SleqpIterate* iterate = solver->iterate;

  assert(sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_DYNAMIC);

  const bool quadratic_model
    = sleqp_settings_bool_value(settings,
                                SLEQP_SETTINGS_BOOL_USE_QUADRATIC_MODEL);

  const bool perform_newton_step
    = quadratic_model
      && sleqp_settings_bool_value(settings,
                                   SLEQP_SETTINGS_BOOL_PERFORM_NEWTON_STEP);

  while (true)
  {
    double current_error_estimate = 0.;

    SLEQP_CALL(sleqp_dyn_func_error_estimate(func, &current_error_estimate));

    const double model_reduction
      = solver->current_merit_value - (*model_trial_value);

    const double required_error_bound
      = compute_required_error_bound(solver, model_reduction);

    if (current_error_estimate <= required_error_bound)
    {
      break;
    }

    sleqp_log_debug("Current accuracy of %e is insufficient, reducing to %e",
                    current_error_estimate,
                    required_error_bound);

    SLEQP_CALL(refine_iterate(solver, problem, iterate, required_error_bound));

    SLEQP_CALL(sleqp_merit_func(solver->merit,
                                iterate,
                                solver->penalty_parameter,
                                &solver->current_merit_value));

    if (perform_newton_step)
    {
      bool failed_eqp_step;
      SLEQP_CALL(compute_trial_point_newton(solver,
                                            trial_iterate,
                                            model_trial_value,
                                            &failed_eqp_step,
                                            full_step));
    }
    else
    {
      SLEQP_CALL(compute_trial_point_simple(solver,
                                            trial_iterate,
                                            model_trial_value,
                                            quadratic_model,
                                            full_step));
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_trial_point_dynamic(SleqpTrialPointSolver* solver,
                            SleqpIterate* trial_iterate,
                            double* trial_merit_value,
                            bool* failed_eqp_step,
                            bool* full_step)
{
  SLEQP_CALL(compute_trial_point_deterministic(solver,
                                               trial_iterate,
                                               trial_merit_value,
                                               failed_eqp_step,
                                               full_step));

  SLEQP_CALL(
    solver_refine_step(solver, trial_iterate, trial_merit_value, full_step));

  SleqpVec* trial_step = sleqp_direction_primal(solver->trial_direction);

  SLEQP_CALL(
    compute_trial_iterate_from_step(solver, trial_step, trial_iterate));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_trial_point_solver_compute_trial_point(SleqpTrialPointSolver* solver,
                                             SleqpIterate* trial_iterate,
                                             double* trial_merit_value,
                                             bool* failed_eqp_step,
                                             bool* full_step,
                                             bool* reject)
{
  assert(solver->trust_radius != SLEQP_NONE);
  assert(solver->lp_trust_radius != SLEQP_NONE);
  assert(solver->penalty_parameter != SLEQP_NONE);

  SleqpProblem* problem = solver->problem;

  SleqpFunc* func = sleqp_problem_func(problem);

  *failed_eqp_step = false;
  // TODO: enable manual rejects
  *reject = false;

  SLEQP_CALL(sleqp_timer_start(solver->elapsed_timer));

  if (sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_DYNAMIC)
  {
    SLEQP_CALL(compute_trial_point_dynamic(solver,
                                           trial_iterate,
                                           trial_merit_value,
                                           failed_eqp_step,
                                           full_step));
  }
  else
  {
    SLEQP_CALL(compute_trial_point_deterministic(solver,
                                                 trial_iterate,
                                                 trial_merit_value,
                                                 failed_eqp_step,
                                                 full_step));
  }

  SLEQP_CALL(sleqp_timer_stop(solver->elapsed_timer));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_trial_point_soc_deterministic(SleqpTrialPointSolver* solver,
                                      SleqpIterate* trial_iterate,
                                      bool* reject)
{
  SleqpIterate* iterate = solver->iterate;

  *reject = false;

  SleqpVec* trial_step = sleqp_direction_primal(solver->trial_direction);
  SleqpVec* soc_step   = sleqp_direction_primal(solver->soc_direction);

  SLEQP_CALL(sleqp_soc_compute_step(solver->soc_data,
                                    solver->aug_jac,
                                    iterate,
                                    trial_step,
                                    trial_iterate,
                                    soc_step));

  SLEQP_CALL(compute_trial_iterate_from_step(solver, soc_step, trial_iterate));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_trial_point_soc_dynamic(SleqpTrialPointSolver* solver,
                                SleqpIterate* trial_iterate,
                                bool* reject)
{
  SleqpProblem* problem = solver->problem;
  SleqpFunc* func       = sleqp_problem_func(problem);

  assert(sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_DYNAMIC);

  SleqpIterate* iterate = solver->iterate;

  SleqpSettings* settings = solver->settings;

  SleqpVec* soc_step = sleqp_direction_primal(solver->soc_direction);

  const double zero_eps
    = sleqp_settings_real_value(solver->settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  const bool quadratic_model
    = sleqp_settings_bool_value(settings,
                                SLEQP_SETTINGS_BOOL_USE_QUADRATIC_MODEL);

  SLEQP_CALL(
    sleqp_soc_compute_step(solver->soc_data,
                           solver->aug_jac,
                           iterate,
                           sleqp_direction_primal(solver->trial_direction),
                           trial_iterate,
                           soc_step));

  SLEQP_CALL(sleqp_direction_reset(solver->soc_direction,
                                   problem,
                                   iterate,
                                   solver->multipliers,
                                   solver->dense_cache,
                                   zero_eps));

  double soc_model_merit = SLEQP_NONE;

  if (quadratic_model)
  {
    SLEQP_CALL(sleqp_merit_quadratic(solver->merit,
                                     iterate,
                                     solver->soc_direction,
                                     solver->penalty_parameter,
                                     &soc_model_merit));
  }
  else
  {
    SLEQP_CALL(sleqp_merit_linear(solver->merit,
                                  iterate,
                                  solver->soc_direction,
                                  solver->penalty_parameter,
                                  &soc_model_merit));
  }

  const double model_reduction = solver->current_merit_value - soc_model_merit;

  double current_error_estimate;

  SLEQP_CALL(sleqp_dyn_func_error_estimate(func, &current_error_estimate));

  const double required_accuracy
    = compute_required_error_bound(solver, model_reduction);

  *reject = (current_error_estimate > required_accuracy);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_trial_point_solver_compute_trial_point_soc(SleqpTrialPointSolver* solver,
                                                 SleqpIterate* trial_iterate,
                                                 bool* reject)
{
  SleqpProblem* problem = solver->problem;
  SleqpFunc* func       = sleqp_problem_func(problem);

  SLEQP_CALL(sleqp_timer_start(solver->elapsed_timer));

  if (sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_DYNAMIC)
  {
    SLEQP_CALL(compute_trial_point_soc_dynamic(solver, trial_iterate, reject));
  }
  else
  {
    SLEQP_CALL(
      compute_trial_point_soc_deterministic(solver, trial_iterate, reject));
  }

  SLEQP_CALL(sleqp_timer_stop(solver->elapsed_timer));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
trial_point_solver_free(SleqpTrialPointSolver** star)
{
  SleqpTrialPointSolver* solver = *star;

  sleqp_free(&solver->cons_weights);

  SLEQP_CALL(sleqp_timer_free(&solver->elapsed_timer));

  sleqp_free(&solver->dense_cache);

  SLEQP_CALL(sleqp_parametric_solver_release(&solver->parametric_solver));

  SLEQP_CALL(
    sleqp_working_set_release(&solver->parametric_original_working_set));

  SLEQP_CALL(sleqp_soc_data_release(&solver->soc_data));

  SLEQP_CALL(sleqp_eqp_solver_release(&solver->eqp_solver));

  SLEQP_CALL(sleqp_working_step_release(&solver->working_step));

  SLEQP_CALL(sleqp_linesearch_release(&solver->linesearch));

  SLEQP_CALL(sleqp_aug_jac_release(&solver->aug_jac));

  SLEQP_CALL(sleqp_fact_release(&solver->fact));

  SLEQP_CALL(sleqp_dual_estimation_release(&solver->estimation_data));

  SLEQP_CALL(sleqp_cauchy_release(&solver->cauchy_data));

  SLEQP_CALL(sleqp_iterate_release(&solver->iterate));

  SLEQP_CALL(sleqp_merit_release(&solver->merit));

  SLEQP_CALL(sleqp_vec_free(&solver->initial_trial_point));
  SLEQP_CALL(sleqp_vec_free(&solver->multipliers));

  SLEQP_CALL(sleqp_direction_release(&solver->trial_direction));
  SLEQP_CALL(sleqp_direction_release(&solver->soc_direction));

  SLEQP_CALL(sleqp_direction_release(&solver->newton_direction));

  SLEQP_CALL(sleqp_vec_free(&solver->estimation_residuals));
  SLEQP_CALL(sleqp_direction_release(&solver->cauchy_direction));
  SLEQP_CALL(sleqp_vec_free(&solver->lp_step));

  SLEQP_CALL(sleqp_settings_release(&solver->settings));

  SLEQP_CALL(sleqp_problem_release(&solver->problem));

  sleqp_free(&solver);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_trial_point_solver_capture(SleqpTrialPointSolver* solver)
{
  ++solver->refcount;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_trial_point_solver_release(SleqpTrialPointSolver** star)
{
  SleqpTrialPointSolver* trial_point_solver = *star;

  if (!trial_point_solver)
  {
    return SLEQP_OKAY;
  }

  if (--trial_point_solver->refcount == 0)
  {
    SLEQP_CALL(trial_point_solver_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
