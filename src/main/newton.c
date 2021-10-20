#include "newton.h"

#include <math.h>
#include <trlib.h>

#include "fail.h"
#include "cmp.h"
#include "feas.h"
#include "iterate.h"
#include "log.h"
#include "mem.h"
#include "timer.h"
#include "util.h"
#include "working_set.h"

#include "sparse/sparse_matrix.h"

#include "tr/tr_solver.h"
#include "tr/trlib_solver.h"
#include "tr/steihaug_solver.h"

struct SleqpNewtonSolver
{
  int refcount;
  SleqpProblem* problem;
  SleqpWorkingStep* working_step;

  SleqpParams* params;
  SleqpOptions* options;

  SleqpIterate* iterate;
  SleqpAugJac* aug_jac;
  double penalty_parameter;

  SleqpSparseVec* gradient;

  SleqpSparseVec* initial_hessian_product;
  SleqpSparseVec* jacobian_product;

  SleqpSparseVec* sparse_cache;
  SleqpSparseVec* tr_step;
  SleqpSparseVec* tr_hessian_product;

  double* dense_cache;

  SleqpTRSolver* tr_solver;

  SleqpTimer* timer;
};

SLEQP_RETCODE sleqp_newton_solver_create(SleqpNewtonSolver** star,
                                         SleqpProblem* problem,
                                         SleqpWorkingStep* step,
                                         SleqpParams* params,
                                         SleqpOptions* options)
{
  SLEQP_CALL(sleqp_malloc(star));

  const int num_variables = sleqp_problem_num_variables(problem);
  const int num_constraints = sleqp_problem_num_constraints(problem);

  SleqpNewtonSolver* solver = *star;

  *solver = (SleqpNewtonSolver) {0};

  solver->refcount = 1;

  solver->problem = problem;
  SLEQP_CALL(sleqp_problem_capture(solver->problem));

  solver->working_step = step;
  SLEQP_CALL(sleqp_working_step_capture(solver->working_step));

  SLEQP_CALL(sleqp_params_capture(params));
  solver->params = params;

  SLEQP_CALL(sleqp_options_capture(options));
  solver->options = options;

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->gradient,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->initial_hessian_product,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->jacobian_product,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->sparse_cache,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->tr_step,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->tr_hessian_product,
                                              num_variables));

  SLEQP_CALL(sleqp_alloc_array(&solver->dense_cache,
                               SLEQP_MAX(num_variables, num_constraints)));

  SLEQP_TR_SOLVER tr_solver = sleqp_options_get_int(options, SLEQP_OPTION_INT_TR_SOLVER);

  if(tr_solver == SLEQP_TR_SOLVER_AUTO)
  {
    SleqpFunc* func = sleqp_problem_func(problem);

    if(sleqp_func_has_psd_hessian(func))
    {
      tr_solver = SLEQP_TR_SOLVER_CG;
    }
    else
    {
      tr_solver = SLEQP_TR_SOLVER_TRLIB;
    }
  }

  if(tr_solver == SLEQP_TR_SOLVER_CG)
  {
    SLEQP_CALL(sleqp_steihaug_solver_create(&solver->tr_solver,
                                            problem,
                                            params,
                                            options));
  }
  else
  {
    // assert(tr_solver == SLEQP_TR_SOLVER_TRLIB);

    SLEQP_CALL(sleqp_trlib_solver_create(&solver->tr_solver,
                                         problem,
                                         params,
                                         options));
  }

  SLEQP_CALL(sleqp_timer_create(&(solver->timer)));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_newton_set_time_limit(SleqpNewtonSolver* solver,
                                          double time_limit)
{
  return sleqp_tr_solver_set_time_limit(solver->tr_solver,
                                        time_limit);
}

SleqpTimer* sleqp_newton_get_timer(SleqpNewtonSolver* solver)
{
  return solver->timer;
}

SLEQP_RETCODE sleqp_newton_set_iterate(SleqpNewtonSolver* solver,
                                       SleqpIterate* iterate,
                                       SleqpAugJac* jacobian,
                                       double trust_radius,
                                       double penalty_parameter)
{
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

SLEQP_RETCODE sleqp_newton_add_violated_multipliers(SleqpNewtonSolver* solver,
                                                    SleqpSparseVec* multipliers)
{
  assert(solver->iterate);

  SleqpSparseVec* cons_dual = sleqp_iterate_get_cons_dual(solver->iterate);

  SleqpSparseVec* violated_cons_mult = sleqp_working_step_get_violated_cons_multipliers(solver->working_step);

  const double zero_eps = sleqp_params_get(solver->params,
                                           SLEQP_PARAM_ZERO_EPS);

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(cons_dual,
                                            violated_cons_mult,
                                            1.,
                                            solver->penalty_parameter,
                                            zero_eps,
                                            multipliers));

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE projection_residuum(SleqpNewtonSolver* solver,
                                  SleqpSparseVec* tr_step,
                                  double* residuum)
{
  const double zero_eps = sleqp_params_get(solver->params,
                                           SLEQP_PARAM_ZERO_EPS);

  SleqpAugJac* jacobian = solver->aug_jac;

  SleqpSparseVec* sparse_cache = solver->sparse_cache;
  SleqpSparseVec* residuals = solver->tr_hessian_product;

  SLEQP_CALL(sleqp_aug_jac_projection(jacobian,
                                      tr_step,
                                      sparse_cache,
                                      NULL));

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(sparse_cache,
                                            tr_step,
                                            1.,
                                            -1.,
                                            zero_eps,
                                            residuals));

  *residuum = sleqp_sparse_vector_inf_norm(residuals);

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE stationarity_residuum(SleqpNewtonSolver* solver,
                                    const SleqpSparseVec* multipliers,
                                    const SleqpSparseVec* gradient,
                                    SleqpSparseVec* tr_step,
                                    double tr_dual,
                                    double* residuum)
{
  SleqpProblem* problem = solver->problem;

  SleqpAugJac* jacobian = solver->aug_jac;
  SleqpSparseVec* sparse_cache = solver->sparse_cache;
  SleqpSparseVec* tr_prod = solver->tr_hessian_product;

  const double zero_eps = sleqp_params_get(solver->params,
                                           SLEQP_PARAM_ZERO_EPS);

  double one = 1.;

  SLEQP_CALL(sleqp_problem_hess_prod(problem,
                                     &one,
                                     tr_step,
                                     multipliers,
                                     tr_prod));

  SLEQP_CALL(sleqp_sparse_vector_add(tr_prod,
                                     gradient,
                                     zero_eps,
                                     sparse_cache));

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(sparse_cache,
                                            tr_step,
                                            1.,
                                            tr_dual,
                                            zero_eps,
                                            tr_prod));

  SLEQP_CALL(sleqp_aug_jac_projection(jacobian,
                                      tr_prod,
                                      sparse_cache,
                                      NULL));

  (*residuum) = sleqp_sparse_vector_inf_norm(sparse_cache);

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE print_residuals(SleqpNewtonSolver* solver,
                              const SleqpSparseVec* multipliers,
                              const SleqpSparseVec* gradient,
                              SleqpSparseVec* tr_step,
                              double trust_radius,
                              double tr_dual)
{
  const double step_norm = sleqp_sparse_vector_norm(tr_step);
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

  sleqp_log_debug("Trust region feasibility residuum: %.14e, stationarity residuum: %.14e",
                  SLEQP_MAX(radius_res, proj_res),
                  stat_res);

  const double stat_tol = sleqp_params_get(solver->params,
                                           SLEQP_PARAM_STATIONARITY_TOL);

  if(stat_res >= stat_tol)
  {
    sleqp_log_warn("Newton stationarity residuum of %e exceeds desired stationarity tolerance of %e",
                   stat_res,
                   stat_tol);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE check_spectrum(SleqpNewtonSolver* solver)
{
  SleqpParams* params = solver->params;

  SleqpProblem* problem = solver->problem;
  SleqpFunc* func = sleqp_problem_func(problem);

  const double eps = sleqp_params_get(params, SLEQP_PARAM_EPS);

  double min_rayleigh = 0., max_rayleigh = 0.;

  SLEQP_CALL(sleqp_newton_current_rayleigh(solver,
                                           &min_rayleigh,
                                           &max_rayleigh));

  assert(min_rayleigh <= max_rayleigh);

  sleqp_log_debug("Spectrum: %.14e, %.14e", min_rayleigh, max_rayleigh);

  if(sleqp_func_has_psd_hessian(func) && sleqp_is_neg(min_rayleigh, eps))
  {
    sleqp_log_warn("Encountered negative Rayleigh quotient (%.14e) on PSD Hessian",
                   min_rayleigh);
  }

  return SLEQP_OKAY;
}

// compute the EQP gradient. Given as the sum of the
// EQP Hessian with the initial solution, the objective
// function gradient and the violated multipliers
static
SLEQP_RETCODE compute_gradient(SleqpNewtonSolver* solver,
                               SleqpSparseVec* multipliers)
{
  assert(solver->iterate);

  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;

  const double zero_eps = sleqp_params_get(solver->params,
                                           SLEQP_PARAM_ZERO_EPS);

  const double penalty_parameter = solver->penalty_parameter;

  SleqpSparseMatrix* cons_jac = sleqp_iterate_get_cons_jac(iterate);
  SleqpSparseVec* func_grad = sleqp_iterate_get_func_grad(iterate);

  double one = 1.;

  SleqpSparseVec* initial_step = sleqp_working_step_get_step(solver->working_step);

  SleqpSparseVec* violated_cons_mult = sleqp_working_step_get_violated_cons_multipliers(solver->working_step);

  SLEQP_CALL(sleqp_problem_hess_prod(problem,
                                     &one,
                                     initial_step,
                                     multipliers,
                                     solver->initial_hessian_product));

  SLEQP_CALL(sleqp_sparse_vector_add(solver->initial_hessian_product,
                                     func_grad,
                                     zero_eps,
                                     solver->sparse_cache));

  SLEQP_CALL(sleqp_sparse_matrix_trans_vector_product(cons_jac,
                                                      violated_cons_mult,
                                                      zero_eps,
                                                      solver->jacobian_product));

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(solver->sparse_cache,
                                            solver->jacobian_product,
                                            1.,
                                            penalty_parameter,
                                            zero_eps,
                                            solver->gradient));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_newton_objective(SleqpNewtonSolver* solver,
                                     const SleqpSparseVec* multipliers,
                                     const SleqpSparseVec* direction,
                                     double* objective)
{
  assert(solver->iterate);

  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;

  const double one = 1.;
  double bilinear_product = 0.;

  SLEQP_CALL(sleqp_problem_hess_bilinear(problem,
                                         &one,
                                         direction,
                                         multipliers,
                                         &bilinear_product));

  SleqpSparseVec* func_grad = sleqp_iterate_get_func_grad(iterate);

  double obj_inner_prod = 0.;

  SLEQP_CALL(sleqp_sparse_vector_dot(direction, func_grad, &obj_inner_prod));

  double cons_inner_prod = 0.;

  SLEQP_CALL(sleqp_sparse_vector_dot(direction,
                                     solver->jacobian_product,
                                     &cons_inner_prod));

  const double offset = sleqp_working_step_get_objective_offset(solver->working_step,
                                                                solver->penalty_parameter);

  *objective = offset + obj_inner_prod + cons_inner_prod + .5*bilinear_product;

  return SLEQP_OKAY;
}

SLEQP_RETCODE print_objective(SleqpNewtonSolver* solver,
                              SleqpSparseVec* multipliers,
                              SleqpSparseVec* newton_step)
{
  double objective;

  SLEQP_CALL(sleqp_newton_objective(solver,
                                    multipliers,
                                    newton_step,
                                    &objective));

  sleqp_log_debug("Newton objective: %g", objective);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_newton_compute_step(SleqpNewtonSolver* solver,
                                        SleqpSparseVec* multipliers,
                                        SleqpSparseVec* newton_step)
{
  assert(solver->iterate);

  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;

  SleqpSparseVec* tr_step = solver->tr_step;

  SleqpAugJac* jacobian = solver->aug_jac;
  double tr_dual = 0.;

  const double eps = sleqp_params_get(solver->params,
                                      SLEQP_PARAM_EPS);

  const double zero_eps = sleqp_params_get(solver->params,
                                           SLEQP_PARAM_ZERO_EPS);

  const double reduced_trust_radius = sleqp_working_step_get_reduced_trust_radius(solver->working_step);

  SleqpSparseVec* initial_step = sleqp_working_step_get_step(solver->working_step);

  // in this case the only feasible solution is the zero vector
  if(sleqp_is_zero(reduced_trust_radius, zero_eps))
  {
    SLEQP_CALL(sleqp_sparse_vector_copy(initial_step, newton_step));

    return SLEQP_OKAY;
  }

  SleqpFunc* func = sleqp_problem_func(problem);

  SleqpTimer* hess_timer = sleqp_func_get_hess_timer(func);

  SleqpTimer* solve_timer = sleqp_aug_jac_solution_timer(solver->aug_jac);

  const double hess_before = sleqp_timer_elapsed(hess_timer);
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

  SLEQP_CALL(sleqp_sparse_vector_add(tr_step,
                                     initial_step,
                                     zero_eps,
                                     newton_step));

  SLEQP_CALL(print_objective(solver, multipliers, newton_step));

  SLEQP_CALL(print_residuals(solver,
                             multipliers,
                             solver->gradient,
                             tr_step,
                             reduced_trust_radius,
                             tr_dual));

#if !defined(NDEBUG)

  // Initial direction and trust region direction
  // must be orthogonal
  {
    double direction_dot;

    SLEQP_CALL(sleqp_sparse_vector_dot(tr_step,
                                       initial_step,
                                       &direction_dot));

    sleqp_assert_is_zero(direction_dot, eps);
  }


  if(sleqp_working_step_in_working_set(solver->working_step))
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

  const double hess_elapsed = sleqp_timer_elapsed(hess_timer) - hess_before;
  const double solve_elapsed = sleqp_timer_elapsed(solve_timer) - solve_before;

  SLEQP_CALL(sleqp_timer_add(solver->timer, -(hess_elapsed + solve_elapsed)));

  return SLEQP_OKAY;
}


SLEQP_RETCODE sleqp_newton_current_rayleigh(SleqpNewtonSolver* solver,
                                            double* min_rayleigh,
                                            double* max_rayleigh)
{
  SLEQP_CALL(sleqp_tr_solver_current_rayleigh(solver->tr_solver,
                                              min_rayleigh,
                                              max_rayleigh));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE newton_data_free(SleqpNewtonSolver** star)
{
  SleqpNewtonSolver* solver = *star;

  if(!solver)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_timer_free(&(solver->timer)));

  SLEQP_CALL(sleqp_tr_solver_release(&solver->tr_solver));

  sleqp_free(&solver->dense_cache);

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->tr_hessian_product));
  SLEQP_CALL(sleqp_sparse_vector_free(&solver->tr_step));
  SLEQP_CALL(sleqp_sparse_vector_free(&solver->sparse_cache));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->jacobian_product));
  SLEQP_CALL(sleqp_sparse_vector_free(&solver->initial_hessian_product));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->gradient));

  SLEQP_CALL(sleqp_aug_jac_release(&solver->aug_jac));
  SLEQP_CALL(sleqp_iterate_release(&solver->iterate));

  SLEQP_CALL(sleqp_options_release(&solver->options));
  SLEQP_CALL(sleqp_params_release(&solver->params));

  SLEQP_CALL(sleqp_working_step_release(&solver->working_step));

  SLEQP_CALL(sleqp_problem_release(&solver->problem));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_newton_solver_capture(SleqpNewtonSolver* data)
{
  ++data->refcount;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_newton_solver_release(SleqpNewtonSolver** star)
{
  SleqpNewtonSolver* data = *star;

  if(!data)
  {
    return SLEQP_OKAY;
  }

  if(--data->refcount == 0)
  {
    SLEQP_CALL(newton_data_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
