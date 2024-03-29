#include "newton.h"

#include <math.h>
#include <trlib.h>

#include "cmp.h"
#include "direction.h"
#include "fail.h"
#include "feas.h"
#include "iterate.h"
#include "log.h"
#include "mem.h"
#include "timer.h"
#include "util.h"
#include "working_set.h"

#include "sparse/mat.h"

#include "tr/steihaug_solver.h"
#include "tr/tr_solver.h"
#include "tr/trlib_solver.h"

typedef struct
{
  int refcount;
  SleqpProblem* problem;
  SleqpWorkingStep* working_step;

  SleqpSettings* settings;

  SleqpIterate* iterate;
  SleqpAugJac* aug_jac;
  double penalty_parameter;

  SleqpVec* gradient;

  SleqpVec* initial_hessian_product;
  SleqpVec* jacobian_product;

  SleqpVec* sparse_cache;
  SleqpVec* tr_step;
  SleqpVec* tr_hessian_product;

  double* dense_cache;

  SleqpTRSolver* tr_solver;

  SleqpTimer* timer;
} NewtonSolver;

static SLEQP_RETCODE
newton_solver_create(NewtonSolver** star,
                     SleqpProblem* problem,
                     SleqpWorkingStep* step,
                     SleqpSettings* settings)
{
  SLEQP_CALL(sleqp_malloc(star));

  const int num_variables   = sleqp_problem_num_vars(problem);
  const int num_constraints = sleqp_problem_num_cons(problem);

  NewtonSolver* solver = *star;

  *solver = (NewtonSolver){0};

  solver->refcount = 1;

  solver->problem = problem;
  SLEQP_CALL(sleqp_problem_capture(solver->problem));

  solver->working_step = step;
  SLEQP_CALL(sleqp_working_step_capture(solver->working_step));

  SLEQP_CALL(sleqp_settings_capture(settings));
  solver->settings = settings;

  SLEQP_CALL(sleqp_vec_create_empty(&solver->gradient, num_variables));

  SLEQP_CALL(
    sleqp_vec_create_empty(&solver->initial_hessian_product, num_variables));

  SLEQP_CALL(sleqp_vec_create_empty(&solver->jacobian_product, num_variables));

  SLEQP_CALL(sleqp_vec_create_empty(&solver->sparse_cache, num_variables));

  SLEQP_CALL(sleqp_vec_create_empty(&solver->tr_step, num_variables));

  SLEQP_CALL(
    sleqp_vec_create_empty(&solver->tr_hessian_product, num_variables));

  SLEQP_CALL(sleqp_alloc_array(&solver->dense_cache,
                               SLEQP_MAX(num_variables, num_constraints)));

  SLEQP_TR_SOLVER tr_solver
    = sleqp_settings_enum_value(settings, SLEQP_SETTINGS_ENUM_TR_SOLVER);

  if (tr_solver == SLEQP_TR_SOLVER_AUTO)
  {
    SleqpFunc* func = sleqp_problem_func(problem);

    if (sleqp_func_has_flags(func, SLEQP_FUNC_HESS_PSD))
    {
      tr_solver = SLEQP_TR_SOLVER_CG;
    }
    else
    {
      tr_solver = SLEQP_TR_SOLVER_TRLIB;
    }
  }

  if (tr_solver == SLEQP_TR_SOLVER_CG)
  {
    SLEQP_CALL(sleqp_steihaug_solver_create(&solver->tr_solver,
                                            problem,
                                            settings));
  }
  else
  {
    // assert(tr_solver == SLEQP_TR_SOLVER_TRLIB);

    SLEQP_CALL(
      sleqp_trlib_solver_create(&solver->tr_solver, problem, settings));
  }

  SLEQP_CALL(sleqp_timer_create(&(solver->timer)));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
newton_solver_current_rayleigh(double* min_rayleigh,
                               double* max_rayleigh,
                               void* data)
{
  NewtonSolver* solver = (NewtonSolver*)data;

  SLEQP_CALL(sleqp_tr_solver_current_rayleigh(solver->tr_solver,
                                              min_rayleigh,
                                              max_rayleigh));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
newton_solver_set_time_limit(double time_limit, void* data)
{
  NewtonSolver* solver = (NewtonSolver*)data;

  return sleqp_tr_solver_set_time_limit(solver->tr_solver, time_limit);
}

SleqpTimer*
sleqp_newton_get_timer(NewtonSolver* solver)
{
  return solver->timer;
}

static SLEQP_RETCODE
newton_solver_set_iterate(SleqpIterate* iterate,
                          SleqpAugJac* jacobian,
                          double trust_radius,
                          double penalty_parameter,
                          void* data)
{
  NewtonSolver* solver = (NewtonSolver*)data;

  solver->penalty_parameter = penalty_parameter;

  {
    SLEQP_CALL(sleqp_iterate_release(&solver->iterate));

    SLEQP_CALL(sleqp_iterate_capture(iterate));

    solver->iterate = iterate;
  }

  {
    SLEQP_CALL(sleqp_aug_jac_release(&solver->aug_jac));

    SLEQP_CALL(sleqp_aug_jac_capture(jacobian));

    solver->aug_jac = jacobian;
  }

  SLEQP_CALL(sleqp_working_step_set_iterate(solver->working_step,
                                            iterate,
                                            jacobian,
                                            trust_radius));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
newton_solver_add_violated_multipliers(SleqpVec* multipliers, void* data)
{
  NewtonSolver* solver = (NewtonSolver*)data;

  assert(solver->iterate);

  SleqpVec* cons_dual = sleqp_iterate_cons_dual(solver->iterate);

  SleqpVec* violated_cons_mult
    = sleqp_working_step_violated_cons_multipliers(solver->working_step);

  const double zero_eps
    = sleqp_settings_real_value(solver->settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  SLEQP_CALL(sleqp_vec_add_scaled(cons_dual,
                                  violated_cons_mult,
                                  1.,
                                  solver->penalty_parameter,
                                  zero_eps,
                                  multipliers));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
projection_residuum(NewtonSolver* solver, SleqpVec* tr_step, double* residuum)
{
  const double zero_eps
    = sleqp_settings_real_value(solver->settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  SleqpAugJac* jacobian = solver->aug_jac;

  SleqpVec* sparse_cache = solver->sparse_cache;
  SleqpVec* residuals    = solver->tr_hessian_product;

  SLEQP_CALL(sleqp_aug_jac_project_nullspace(jacobian, tr_step, sparse_cache));

  SLEQP_CALL(
    sleqp_vec_add_scaled(sparse_cache, tr_step, 1., -1., zero_eps, residuals));

  *residuum = sleqp_vec_inf_norm(residuals);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
stationarity_residuum(NewtonSolver* solver,
                      const SleqpVec* multipliers,
                      const SleqpVec* gradient,
                      SleqpVec* tr_step,
                      double tr_dual,
                      double* residuum)
{
  SleqpProblem* problem = solver->problem;

  SleqpAugJac* jacobian  = solver->aug_jac;
  SleqpVec* sparse_cache = solver->sparse_cache;
  SleqpVec* tr_prod      = solver->tr_hessian_product;

  const double zero_eps
    = sleqp_settings_real_value(solver->settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  SLEQP_CALL(
    sleqp_problem_hess_prod(problem, tr_step, multipliers, tr_prod));

  SLEQP_CALL(sleqp_vec_add(tr_prod, gradient, zero_eps, sparse_cache));

  SLEQP_CALL(sleqp_vec_add_scaled(sparse_cache,
                                  tr_step,
                                  1.,
                                  tr_dual,
                                  zero_eps,
                                  tr_prod));

  SLEQP_CALL(sleqp_aug_jac_project_nullspace(jacobian, tr_prod, sparse_cache));

  (*residuum) = sleqp_vec_inf_norm(sparse_cache);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
print_residuals(NewtonSolver* solver,
                const SleqpVec* multipliers,
                const SleqpVec* gradient,
                SleqpVec* tr_step,
                double trust_radius,
                double tr_dual)
{
  const double step_norm  = sleqp_vec_norm(tr_step);
  const double radius_res = SLEQP_MAX(step_norm - trust_radius, 0.);

  double proj_res = 0.;

  SLEQP_CALL(projection_residuum(solver, tr_step, &proj_res));

  double stat_res = 0.;

  SLEQP_CALL(stationarity_residuum(solver,
                                   multipliers,
                                   gradient,
                                   tr_step,
                                   tr_dual,
                                   &stat_res));

  sleqp_log_debug(
    "Trust region feasibility residuum: %.14e, stationarity residuum: %.14e",
    SLEQP_MAX(radius_res, proj_res),
    stat_res);

  const double stat_tol
    = sleqp_settings_real_value(solver->settings, SLEQP_SETTINGS_REAL_STAT_TOL);

  if (stat_res >= stat_tol)
  {
    sleqp_log_warn("Newton stationarity residuum of %e exceeds desired "
                   "stationarity tolerance of %e",
                   stat_res,
                   stat_tol);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
check_spectrum(NewtonSolver* solver)
{
  SleqpSettings* settings = solver->settings;

  SleqpProblem* problem = solver->problem;
  SleqpFunc* func       = sleqp_problem_func(problem);

  const double eps = sleqp_settings_real_value(settings, SLEQP_SETTINGS_REAL_EPS);

  double min_rayleigh = 0., max_rayleigh = 0.;

  SLEQP_CALL(
    newton_solver_current_rayleigh(&min_rayleigh, &max_rayleigh, solver));

  assert(min_rayleigh <= max_rayleigh);

  sleqp_log_debug("Spectrum: %.14e, %.14e", min_rayleigh, max_rayleigh);

  if (sleqp_func_has_flags(func, SLEQP_FUNC_HESS_PSD)
      && sleqp_is_neg(min_rayleigh, eps))
  {
    sleqp_log_warn(
      "Encountered negative Rayleigh quotient (%.14e) on PSD Hessian",
      min_rayleigh);
  }

  return SLEQP_OKAY;
}

// compute the EQP gradient. Given as the sum of the
// EQP Hessian with the initial solution, the objective
// function gradient and the violated multipliers
static SLEQP_RETCODE
compute_gradient(NewtonSolver* solver, const SleqpVec* multipliers)
{
  assert(solver->iterate);

  SleqpIterate* iterate = solver->iterate;

  SleqpDirection* initial_direction
    = sleqp_working_step_direction(solver->working_step);

  const double zero_eps
    = sleqp_settings_real_value(solver->settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  const double penalty_parameter = solver->penalty_parameter;

  SleqpMat* cons_jac = sleqp_iterate_cons_jac(iterate);
  SleqpVec* obj_grad = sleqp_iterate_obj_grad(iterate);

  SleqpVec* violated_cons_mult
    = sleqp_working_step_violated_cons_multipliers(solver->working_step);

  SleqpVec* initial_hess = sleqp_direction_hess(initial_direction);

  SLEQP_CALL(
    sleqp_vec_add(initial_hess, obj_grad, zero_eps, solver->sparse_cache));

  SLEQP_CALL(sleqp_mat_mult_vec_trans(cons_jac,
                                      violated_cons_mult,
                                      zero_eps,
                                      solver->jacobian_product));

  SLEQP_CALL(sleqp_vec_add_scaled(solver->sparse_cache,
                                  solver->jacobian_product,
                                  1.,
                                  penalty_parameter,
                                  zero_eps,
                                  solver->gradient));

  return SLEQP_OKAY;
}

// Compute Newton objective, assuming that
// the direction produces a linear violation
// equal to that of the initial step into
// the working set Jacobian
static SLEQP_RETCODE
newton_objective(NewtonSolver* solver,
                 const SleqpDirection* direction,
                 double* objective)
{
  assert(solver->iterate);

  *objective = sleqp_working_step_newton_obj_offset(solver->working_step,
                                                    solver->penalty_parameter);

  // add objective gradient
  {
    *objective += *sleqp_direction_obj_grad(direction);
  }

  // add hessian
  {
    SleqpVec* primal = sleqp_direction_primal(direction);
    SleqpVec* hess   = sleqp_direction_hess(direction);

    double bilinear_product;

    SLEQP_CALL(sleqp_vec_dot(primal, hess, &bilinear_product));

    *objective += .5 * bilinear_product;
  }

  assert(sleqp_is_finite(*objective));

  return SLEQP_OKAY;
}

/*
static SLEQP_RETCODE
print_objective(NewtonSolver* solver, const SleqpDirection* direction)
{
  double objective;

  SLEQP_CALL(newton_objective(solver, direction, &objective));

  sleqp_log_debug("Newton objective: %g", objective);

  return SLEQP_OKAY;
}
*/

static SLEQP_RETCODE
newton_solver_compute_direction(const SleqpVec* multipliers,
                                SleqpDirection* newton_direction,
                                void* data)
{
  NewtonSolver* solver = (NewtonSolver*)data;
  assert(solver->iterate);

  SleqpProblem* problem          = solver->problem;
  SleqpIterate* iterate          = solver->iterate;
  SleqpWorkingStep* working_step = solver->working_step;

  SleqpVec* tr_step = solver->tr_step;

  SleqpAugJac* jacobian = solver->aug_jac;
  double tr_dual        = 0.;

  SleqpVec* newton_step = sleqp_direction_primal(newton_direction);

  const double eps = sleqp_settings_real_value(solver->settings, SLEQP_SETTINGS_REAL_EPS);

  const double zero_eps
    = sleqp_settings_real_value(solver->settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  const double reduced_trust_radius
    = sleqp_working_step_reduced_trust_radius(solver->working_step);

  SLEQP_CALL(sleqp_working_step_set_multipliers(working_step, multipliers));

  SleqpDirection* initial_direction
    = sleqp_working_step_direction(solver->working_step);

  {
    double newton_obj;

    SLEQP_CALL(newton_objective(solver, initial_direction, &newton_obj));

    sleqp_log_debug("Newton objective of initial step: %g", newton_obj);
  }

  SleqpVec* initial_step = sleqp_direction_primal(initial_direction);

#if SLEQP_DEBUG

  {
    bool valid;

    SLEQP_CALL(sleqp_direction_check(sleqp_working_step_direction(working_step),
                                     problem,
                                     iterate,
                                     multipliers,
                                     solver->dense_cache,
                                     zero_eps,
                                     &valid));

    sleqp_num_assert(valid);
  }

#endif

  // in this case the only feasible solution is the zero vector
  if (sleqp_is_zero(reduced_trust_radius, zero_eps))
  {
    SLEQP_CALL(sleqp_direction_copy(initial_direction, newton_direction));

    return SLEQP_OKAY;
  }

  SleqpFunc* func = sleqp_problem_func(problem);

  SleqpTimer* hess_timer = sleqp_func_get_hess_timer(func);

  SleqpTimer* solve_timer = sleqp_aug_jac_solution_timer(solver->aug_jac);

  const double hess_before  = sleqp_timer_elapsed(hess_timer);
  const double solve_before = sleqp_timer_elapsed(solve_timer);

  SLEQP_CALL(sleqp_timer_start(solver->timer));

  SLEQP_CALL(compute_gradient(solver, multipliers));

  SLEQP_CALL(sleqp_tr_solver_solve(solver->tr_solver,
                                   jacobian,
                                   multipliers,
                                   solver->gradient,
                                   tr_step,
                                   reduced_trust_radius,
                                   &tr_dual));

  SLEQP_CALL(check_spectrum(solver));

  SLEQP_CALL(sleqp_vec_add(tr_step, initial_step, zero_eps, newton_step));

  SLEQP_CALL(sleqp_direction_reset(newton_direction,
                                   solver->problem,
                                   solver->iterate,
                                   multipliers,
                                   solver->dense_cache,
                                   zero_eps));

  {
    double newton_obj;

    SLEQP_CALL(newton_objective(solver, newton_direction, &newton_obj));

    sleqp_log_debug("Newton objective of Newton step: %g", newton_obj);
  }

  SLEQP_CALL(print_residuals(solver,
                             multipliers,
                             solver->gradient,
                             tr_step,
                             reduced_trust_radius,
                             tr_dual));

#if SLEQP_DEBUG

  // Initial direction and trust region direction
  // must be orthogonal
  {
    double direction_dot;

    SLEQP_CALL(sleqp_vec_dot(tr_step, initial_step, &direction_dot));

    sleqp_assert_is_zero(direction_dot, eps);
  }

  if (sleqp_working_step_in_working_set(solver->working_step))
  {
    // Direction must be in working set
    bool in_working_set = false;

    SLEQP_CALL(sleqp_direction_in_working_set(problem,
                                              iterate,
                                              newton_step,
                                              solver->dense_cache,
                                              eps,
                                              &in_working_set));

    sleqp_num_assert(in_working_set);
  }

#endif

  SLEQP_CALL(sleqp_timer_stop(solver->timer));

  const double hess_elapsed  = sleqp_timer_elapsed(hess_timer) - hess_before;
  const double solve_elapsed = sleqp_timer_elapsed(solve_timer) - solve_before;

  SLEQP_CALL(sleqp_timer_add(solver->timer, -(hess_elapsed + solve_elapsed)));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
newton_solver_free(void* data)
{
  NewtonSolver* solver = (NewtonSolver*)data;

  if (!solver)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_timer_free(&(solver->timer)));

  SLEQP_CALL(sleqp_tr_solver_release(&solver->tr_solver));

  sleqp_free(&solver->dense_cache);

  SLEQP_CALL(sleqp_vec_free(&solver->tr_hessian_product));
  SLEQP_CALL(sleqp_vec_free(&solver->tr_step));
  SLEQP_CALL(sleqp_vec_free(&solver->sparse_cache));

  SLEQP_CALL(sleqp_vec_free(&solver->jacobian_product));
  SLEQP_CALL(sleqp_vec_free(&solver->initial_hessian_product));

  SLEQP_CALL(sleqp_vec_free(&solver->gradient));

  SLEQP_CALL(sleqp_aug_jac_release(&solver->aug_jac));
  SLEQP_CALL(sleqp_iterate_release(&solver->iterate));

  SLEQP_CALL(sleqp_settings_release(&solver->settings));

  SLEQP_CALL(sleqp_working_step_release(&solver->working_step));

  SLEQP_CALL(sleqp_problem_release(&solver->problem));

  sleqp_free(&solver);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_newton_solver_create(SleqpEQPSolver** star,
                           SleqpProblem* problem,
                           SleqpSettings* settings,
                           SleqpWorkingStep* step)
{
  SleqpEQPCallbacks callbacks
    = {.set_iterate              = newton_solver_set_iterate,
       .set_time_limit           = newton_solver_set_time_limit,
       .add_violated_multipliers = newton_solver_add_violated_multipliers,
       .compute_direction        = newton_solver_compute_direction,
       .current_rayleigh         = newton_solver_current_rayleigh,
       .free                     = newton_solver_free};

  NewtonSolver* solver;

  SLEQP_CALL(newton_solver_create(&solver, problem, step, settings));

  SLEQP_CALL(sleqp_eqp_solver_create(star, &callbacks, (void*)solver));

  return SLEQP_OKAY;
}
