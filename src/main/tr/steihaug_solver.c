/*
 * An iterative solver for the trust-region subproblem using projected Conjugate
 * Gradients (CG) with Steihaug's modification for the boundary case. The
 * augmented Jacobian system is used to project onto the nullspace of the active
 * set identified in the LP step. The (1,1)-block of the projector is currently
 * the identity, but could contain a Hessian preconditioner.
 */

#include "steihaug_solver.h"

#include <math.h>

#include "aug_jac/aug_jac.h"
#include "cmp.h"
#include "fail.h"
#include "log.h"
#include "mem.h"
#include "params.h"
#include "sparse/pub_vec.h"
#include "tr/tr_util.h"

static const double tolerance_factor = 1e-2;

typedef struct
{
  SleqpProblem* problem;
  SleqpParams* params;

  double time_limit;
  int max_iter;

  SleqpVec* d;  // projected descent direction
  SleqpVec* Bd; // Hessian-times-direction product
  SleqpVec* g;  // projected residual
  SleqpVec* r;  // full-space residual
  SleqpVec* z;  // CG iterate

  SleqpVec* sparse_cache;

  double min_rayleigh, max_rayleigh;

  SleqpTimer* timer;
} SleqpSteihaugSolver;

static SLEQP_RETCODE
steihaug_solver_free(void** star)
{
  SleqpSteihaugSolver* solver = (SleqpSteihaugSolver*)(*star);

  if (!solver)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_timer_free(&solver->timer));

  SLEQP_CALL(sleqp_vec_free(&solver->sparse_cache));
  SLEQP_CALL(sleqp_vec_free(&solver->z));
  SLEQP_CALL(sleqp_vec_free(&solver->r));
  SLEQP_CALL(sleqp_vec_free(&solver->g));
  SLEQP_CALL(sleqp_vec_free(&solver->Bd));
  SLEQP_CALL(sleqp_vec_free(&solver->d));

  SLEQP_CALL(sleqp_params_release(&solver->params));
  SLEQP_CALL(sleqp_problem_release(&solver->problem));

  sleqp_free(star);

  return SLEQP_OKAY;
}

#if !defined(NDEBUG)

static SLEQP_RETCODE
check_projection(SleqpSteihaugSolver* solver,
                 SleqpAugJac* jacobian,
                 SleqpVec* step)
{
  SleqpVec* sparse_cache = solver->sparse_cache;

  assert(step != sparse_cache);

  const double eps = sleqp_params_value(solver->params, SLEQP_PARAM_EPS);

  SLEQP_NUM_ASSERT_PARAM(eps);

  SLEQP_CALL(sleqp_aug_jac_project_nullspace(jacobian, step, sparse_cache));

  sleqp_num_assert(sleqp_vec_eq(step, sparse_cache, eps));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
check_trust_radius(SleqpSteihaugSolver* solver,
                   double trust_radius,
                   SleqpVec* step)
{
  SleqpVec* sparse_cache = solver->sparse_cache;

  assert(step != sparse_cache);

  const double eps = sleqp_params_value(solver->params, SLEQP_PARAM_EPS);

  SLEQP_NUM_ASSERT_PARAM(eps);

  sleqp_assert_is_leq(sleqp_vec_norm(step), trust_radius, eps);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
check_direction(SleqpSteihaugSolver* solver,
                SleqpAugJac* jacobian,
                double trust_radius,
                SleqpVec* step)
{
  SLEQP_CALL(check_projection(solver, jacobian, step));
  SLEQP_CALL(check_trust_radius(solver, trust_radius, step));

  return SLEQP_OKAY;
}

#endif

static SLEQP_RETCODE
steihaug_collect_rayleigh(SleqpSteihaugSolver* solver,
                          const SleqpVec* direction,
                          const SleqpVec* product)
{
  const double dir_normsq = sleqp_vec_norm_sq(direction);

  if (dir_normsq == 0.)
  {
    return SLEQP_OKAY;
  }

  double dot_product;

  SLEQP_CALL(sleqp_vec_dot(direction, product, &dot_product));

  const double cur_rayleigh = dot_product / dir_normsq;

  solver->min_rayleigh = SLEQP_MIN(solver->min_rayleigh, cur_rayleigh);

  solver->max_rayleigh = SLEQP_MAX(solver->max_rayleigh, cur_rayleigh);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
steihaug_solver_rayleigh(double* min_rayleigh,
                         double* max_rayleigh,
                         void* solver_data)
{
  SleqpSteihaugSolver* solver = (SleqpSteihaugSolver*)solver_data;

  (*min_rayleigh) = solver->min_rayleigh;
  (*max_rayleigh) = solver->max_rayleigh;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
steihaug_solver_solve(SleqpAugJac* jacobian,
                      const SleqpVec* multipliers,
                      const SleqpVec* gradient,
                      SleqpVec* newton_step,
                      double trust_radius,
                      double* tr_dual,
                      double time_limit,
                      void* solver_data)
{
  SleqpSteihaugSolver* solver = (SleqpSteihaugSolver*)solver_data;

  solver->min_rayleigh = 1.;
  solver->max_rayleigh = 1.;

  const double one = 1.;

  const double stat_eps
    = sleqp_params_value(solver->params, SLEQP_PARAM_STAT_TOL);

  const double rel_tol = stat_eps * tolerance_factor;

  SleqpProblem* problem = solver->problem;
  SleqpParams* params   = solver->params;

  *tr_dual = SLEQP_NONE;

  const double eps      = sleqp_params_value(params, SLEQP_PARAM_EPS);
  const double zero_eps = sleqp_params_value(params, SLEQP_PARAM_ZERO_EPS);

  SLEQP_NUM_ASSERT_PARAM(eps);

  const double rel_tol_sq = rel_tol * rel_tol;

  double dBd;
  double alpha;
  double beta;
  double d_nrm_sq;
  double gd;
  double zBd;

  double z_curr_nrm_sq = 0.;
  double z_next_nrm_sq = 0.;

  SLEQP_CALL(sleqp_timer_start(solver->timer));

  SLEQP_CALL(sleqp_vec_clear(newton_step));

  // set z0 such that P[z0] = 0
  SLEQP_CALL(sleqp_vec_clear(solver->z));

  // set r0 = nabla f_k
  SLEQP_CALL(sleqp_vec_copy(gradient, solver->r));

  // set g0 = P[r0]
  SLEQP_CALL(sleqp_aug_jac_project_nullspace(jacobian, solver->r, solver->g));

  // set d0 = -P[r0] = - P[nabla f_k] = -g0
  SLEQP_CALL(sleqp_vec_copy(solver->g, solver->d));
  SLEQP_CALL(sleqp_vec_scale(solver->d, -1.));

  // if ||d0|| < eps_k: return p_k = P[z_0] = 0
  d_nrm_sq = sleqp_vec_norm_sq(solver->d);

  if (d_nrm_sq < rel_tol_sq)
  {
    SLEQP_CALL(sleqp_vec_copy(solver->z, newton_step));
    SLEQP_CALL(sleqp_timer_stop(solver->timer));
    return SLEQP_OKAY;
  }

  // compute r_j^T * g_j
  double r_dot_g;

  SLEQP_CALL(sleqp_vec_dot(solver->r, solver->g, &r_dot_g));

  bool reached_time_limit = false;

  // loop over pCG iterations j
  for (int iteration = 0;; ++iteration)
  {
    if (solver->max_iter != SLEQP_NONE && iteration >= solver->max_iter)
    {
      break;
    }

    if (time_limit != SLEQP_NONE
        && sleqp_timer_elapsed(solver->timer) >= time_limit)
    {
      reached_time_limit = true;
      break;
    }

#if !defined(NDEBUG)
    SLEQP_CALL(check_projection(solver, jacobian, solver->d));
#endif

    // if |r_{j+1}^T * g_{j+1}| < eps_k:
    if (fabs(r_dot_g) < rel_tol_sq)
    {
      sleqp_log_debug("CG solver found interior solution after %d iterations",
                      iteration);

      // return p_k = z_{j+1}
      SLEQP_CALL(sleqp_vec_copy(solver->z, newton_step));

      break;
    }

    // compute B_k * d_j
    SLEQP_CALL(sleqp_problem_hess_prod(problem,
                                       &one,
                                       solver->d,
                                       multipliers,
                                       solver->Bd));

    SLEQP_CALL(steihaug_collect_rayleigh(solver, solver->d, solver->Bd));

    // compute d_j^T * (B_k * d_j)
    SLEQP_CALL(sleqp_vec_dot(solver->d, solver->Bd, &dBd));

    // if d_j^T * B_k * d_j <= 0:
    if (dBd <= 0.)
    {
      // z is a feasible iterate and d is a null-space direction, so P^T[p_k] =
      // 0 find tau such that p_k = z_j + tau * d_j solves min_{tau}   m_k(p) =
      // f_k + nabla f_k^T * p + 0.5 * p^T * B_k * p subject to  ||p_k|| = Delta

      // solving for tau results in tau = - z^T * d / ||d||^2 \pm Delta / ||d||
      double z_dot_d;
      SLEQP_CALL(sleqp_vec_dot(solver->z, solver->d, &z_dot_d));

      const double d_nrm_sq = sleqp_vec_norm_sq(solver->d);
      const double z_nrm_sq = z_curr_nrm_sq;

      assert(d_nrm_sq > 0.);

      const double inner
        = z_dot_d * z_dot_d
          - d_nrm_sq * (z_nrm_sq - trust_radius * trust_radius);

      const double tau_min = 1. / d_nrm_sq * (-z_dot_d - sqrt(inner));
      const double tau_max = 1. / d_nrm_sq * (-z_dot_d + sqrt(inner));

      // compute dot products g^T * d and z^T * (B * d)
      SLEQP_CALL(sleqp_vec_dot(gradient, solver->d, &gd));
      SLEQP_CALL(sleqp_vec_dot(solver->z, solver->Bd, &zBd));

      const double tau_min_obj = tau_min * ((gd + zBd) + 0.5 * tau_min * dBd);
      const double tau_max_obj = tau_max * ((gd + zBd) + 0.5 * tau_max * dBd);

      // pick the boundary intersection with smaller quadratic model value
      // q(p) = tau * (g^T * d + z^T * B * d) + 0.5 * tau^2 * d^T * B * d +
      // const.
      const double tau = (tau_min_obj < tau_max_obj) ? tau_min : tau_max;

      // return p_k
      SLEQP_CALL(sleqp_vec_add_scaled(solver->z,
                                      solver->d,
                                      1.,
                                      tau,
                                      zero_eps,
                                      newton_step));

      sleqp_num_assert(
        sleqp_is_eq(sleqp_vec_norm(newton_step), trust_radius, eps));

      sleqp_log_debug(
        "CG solver found negative curvature direction after %d iterations",
        iteration);

      break;
    }

    // set alpha_j = (r_j^T * g_j) / (d_j^T * B_k * d_j)
    alpha = r_dot_g / dBd;

    // set z_{j+1} = z_j + alpha_j * d_j
    SLEQP_CALL(sleqp_vec_add_scaled(solver->z,
                                    solver->d,
                                    1.,
                                    alpha,
                                    zero_eps,
                                    solver->sparse_cache));

    z_next_nrm_sq = sleqp_vec_norm_sq(solver->sparse_cache);

    // if ||z_{j+1}|| >= Delta_k:

    if (z_next_nrm_sq >= trust_radius * trust_radius)
    {
      SLEQP_CALL(sleqp_tr_compute_bdry_sol(solver->z,
                                           solver->d,
                                           params,
                                           trust_radius,
                                           newton_step));

      sleqp_log_debug("CG solver found boundary solution after %d iterations",
                      iteration);

#if !defined(NDEBUG)
      SLEQP_CALL(check_direction(solver, jacobian, trust_radius, newton_step));
#endif

      break;
    }

    SLEQP_CALL(sleqp_vec_copy(solver->sparse_cache, solver->z));

    z_curr_nrm_sq = z_next_nrm_sq;

    // set r_{j+1} = r_j + alpha_j * (B_k * d_j)
    SLEQP_CALL(sleqp_vec_add_scaled(solver->r,
                                    solver->Bd,
                                    1.,
                                    alpha,
                                    zero_eps,
                                    solver->sparse_cache));

    SLEQP_CALL(sleqp_vec_copy(solver->sparse_cache, solver->r));

    // set g_{j+1} = P[r_{j+1}]
    SLEQP_CALL(sleqp_aug_jac_project_nullspace(jacobian, solver->r, solver->g));

#if !defined(NDEBUG)
    SLEQP_CALL(check_projection(solver, jacobian, solver->g));
#endif

    // set beta_{j+1} = r_{j+1}^T * g_{j+1} / r_j^T * g_j
    beta = 1. / r_dot_g;
    SLEQP_CALL(sleqp_vec_dot(solver->r, solver->g, &r_dot_g));
    beta *= r_dot_g;

    // set d_{j+1} = - g_{j+1} + beta_{j+1} * d_j
    SLEQP_CALL(sleqp_vec_add_scaled(solver->g,
                                    solver->d,
                                    -1.,
                                    beta,
                                    zero_eps,
                                    solver->sparse_cache));

    SLEQP_CALL(sleqp_vec_copy(solver->sparse_cache, solver->d));

#if !defined(NDEBUG)
    SLEQP_CALL(check_projection(solver, jacobian, solver->g));
    SLEQP_CALL(check_projection(solver, jacobian, solver->d));
#endif
  }
  // end loop

  sleqp_num_assert(
    sleqp_is_leq(sleqp_vec_norm(newton_step), trust_radius, eps));

  SLEQP_CALL(sleqp_timer_stop(solver->timer));

  if (reached_time_limit)
  {
    return SLEQP_ABORT_TIME;
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_steihaug_solver_create(SleqpTRSolver** solver_star,
                             SleqpProblem* problem,
                             SleqpParams* params,
                             SleqpOptions* options)
{
  SleqpSteihaugSolver* solver = NULL;

  const int num_variables = sleqp_problem_num_vars(problem);

  SLEQP_CALL(sleqp_malloc(&solver));

  *solver = (SleqpSteihaugSolver){0};

  solver->problem = problem;
  SLEQP_CALL(sleqp_problem_capture(solver->problem));

  SLEQP_CALL(sleqp_params_capture(params));
  solver->params = params;

  solver->max_iter
    = sleqp_options_int_value(options, SLEQP_OPTION_INT_MAX_NEWTON_ITERATIONS);

  SLEQP_CALL(sleqp_vec_create_empty(&solver->d, num_variables));
  SLEQP_CALL(sleqp_vec_create_empty(&solver->Bd, num_variables));
  SLEQP_CALL(sleqp_vec_create_empty(&solver->g, num_variables));
  SLEQP_CALL(sleqp_vec_create_empty(&solver->r, num_variables));
  SLEQP_CALL(sleqp_vec_create_empty(&solver->z, num_variables));
  SLEQP_CALL(sleqp_vec_create_empty(&solver->sparse_cache, num_variables));

  SLEQP_CALL(sleqp_timer_create(&solver->timer));

  SleqpTRCallbacks callbacks = {.solve    = steihaug_solver_solve,
                                .rayleigh = steihaug_solver_rayleigh,
                                .free     = steihaug_solver_free};

  SLEQP_CALL(sleqp_tr_solver_create(solver_star, &callbacks, (void*)solver));

  return SLEQP_OKAY;
}
