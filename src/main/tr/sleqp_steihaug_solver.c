/*
 * An iterative solver for the trust-region subproblem using projected Conjugate Gradients (CG)
 * with Steihaug's modification for the boundary case. The augmented Jacobian system is used to
 * project onto the nullspace of the active set identified in the LP step. The (1,1)-block of the
 * projector is currently the identity, but could contain a Hessian preconditioner.
 */

#include "sleqp_steihaug_solver.h"

#include <math.h>

#include "sleqp_assert.h"
#include "sleqp_cmp.h"
#include "sleqp_mem.h"

typedef struct SleqpSteihaugSolver SleqpSteihaugSolver;


struct SleqpSteihaugSolver
{
  int refcount;
  SleqpProblem* problem;
  SleqpParams* params;

  double time_limit;
  int max_iter;

  SleqpSparseVec* d;            // projected descent direction
  SleqpSparseVec* Bd;           // Hessian-times-direction product
  SleqpSparseVec* g;            // projected residual
  SleqpSparseVec* r;            // full-space residual
  SleqpSparseVec* z;            // CG iterate

  SleqpSparseVec* sparse_cache; // a temporary t for computing t <- a*x + b*y; x <- t

  SleqpTimer* timer;
};


static SLEQP_RETCODE steihaug_solver_free(void **star)
{
  SleqpSteihaugSolver* solver = (SleqpSteihaugSolver*) (*star);

  if (!solver)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_timer_free(&solver->timer));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->sparse_cache));
  SLEQP_CALL(sleqp_sparse_vector_free(&solver->z));
  SLEQP_CALL(sleqp_sparse_vector_free(&solver->r));
  SLEQP_CALL(sleqp_sparse_vector_free(&solver->g));
  SLEQP_CALL(sleqp_sparse_vector_free(&solver->Bd));
  SLEQP_CALL(sleqp_sparse_vector_free(&solver->d));

  SLEQP_CALL(sleqp_params_release (&solver->params));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE steihaug_solver_solve(SleqpAugJacobian* jacobian,
                                    SleqpSparseVec* multipliers,
                                    SleqpSparseVec* gradient,
                                    SleqpSparseVec* newton_step,
                                    double trust_radius,
                                    double time_limit,
                                    void* solver_data)
{
  SleqpSteihaugSolver* solver = (SleqpSteihaugSolver*) solver_data;

  const double one = 1.;
  const double rel_tol = sleqp_params_get (solver->params, SLEQP_PARAM_NEWTON_RELATIVE_TOL);

  SleqpProblem* problem = solver->problem;
  SleqpFunc* func = problem->func;
  SleqpParams* params = solver->params;

  const double eps = sleqp_params_get(params, SLEQP_PARAM_EPS);

  const double rel_tol_sq = rel_tol*rel_tol;

  double dBd;
  double alpha;
  double beta;
  double d_nrm_sq;
  double tau;
  double tau1;
  double tau2;
  double gd;
  double zBd;

  double z_curr_nrm_sq = 0.;
  double z_next_nrm_sq = 0.;

  SLEQP_CALL(sleqp_timer_start(solver->timer));

  SLEQP_CALL(sleqp_sparse_vector_clear(newton_step));

  // set z0 such that P[z0] = 0
  SLEQP_CALL(sleqp_sparse_vector_clear(solver->z));

  // set r0 = nabla f_k
  SLEQP_CALL(sleqp_sparse_vector_copy(gradient, solver->r));

  // set g0 = P[r0]
  SLEQP_CALL(sleqp_aug_jacobian_projection(jacobian, solver->r, solver->g, NULL));

  // set d0 = -P[r0] = - P[nabla f_k] = -g0
  SLEQP_CALL(sleqp_sparse_vector_copy(solver->g, solver->d));
  SLEQP_CALL(sleqp_sparse_vector_scale(solver->d, -1.));

  // if ||d0|| < eps_k: return p_k = P[z_0] = 0
  d_nrm_sq = sleqp_sparse_vector_norm_sq(solver->d);

  if(d_nrm_sq < rel_tol_sq)
  {
    SLEQP_CALL(sleqp_sparse_vector_copy(solver->z, newton_step));
    SLEQP_CALL(sleqp_timer_stop(solver->timer));
    return SLEQP_OKAY;
  }

  // compute r_j^T * g_j
  double r_dot_g;

  SLEQP_CALL(sleqp_sparse_vector_dot(solver->r, solver->g, &r_dot_g));

  // loop over pCG iterations j
  for(int iteration = 0;; ++iteration)
  {
    if(solver->max_iter != SLEQP_NONE && iteration >= solver->max_iter)
    {
      break;
    }

    if(time_limit != SLEQP_NONE && sleqp_timer_elapsed(solver->timer) >= time_limit)
    {
      break;
    }

    // if |r_{j+1}^T * g_{j+1}| < eps_k:
    if(fabs(r_dot_g) < rel_tol_sq)
    {
      sleqp_log_debug("CG solver found interior solution after %d iterations",
                      iteration);

      // return p_k = z_{j+1}
      SLEQP_CALL(sleqp_sparse_vector_copy(solver->z, newton_step));
      break;
    }

    // compute B_k * d_j
    SLEQP_CALL(sleqp_func_hess_prod(func, &one, solver->d, multipliers, solver->Bd));

    // compute d_j^T * (B_k * d_j)
    SLEQP_CALL(sleqp_sparse_vector_dot(solver->d, solver->Bd, &dBd));

    // if d_j^T * B_k * d_j <= 0:
    if (dBd <= 0.)
    {
      // z is a feasible iterate and d is a null-space direction, so P^T[p_k] = 0
      // find tau such that p_k = z_j + tau * d_j solves
      // min_{tau}   m_k(p) = f_k + nabla f_k^T * p + 0.5 * p^T * B_k * p
      // subject to  ||p_k|| = Delta

      // solving for tau results in tau = - z^T * d / ||d||^2 \pm Delta / ||d||
      SLEQP_CALL(sleqp_sparse_vector_dot(solver->z, solver->d, &tau));
      const double d_nrm_sq = sleqp_sparse_vector_norm_sq(solver->d);
      tau1 = - tau / d_nrm_sq + trust_radius / sqrt (d_nrm_sq);        /// @todo verify again
      tau2 = - tau / d_nrm_sq - trust_radius / sqrt (d_nrm_sq);

      // compute dot products g^T * d and z^T * (B * d)
      SLEQP_CALL(sleqp_sparse_vector_dot(gradient, solver->d, &gd));
      SLEQP_CALL(sleqp_sparse_vector_dot(solver->z, solver->Bd, &zBd));

      // pick the boundary intersection with smaller quadratic model value
      // q(p) = tau * (g^T * d + z^T * B * d) + 0.5 * tau^2 * d^T * B * d + const.
      tau = (tau1*((gd + zBd) + 0.5*tau1*dBd) < tau2*((gd + zBd) + 0.5*tau2*dBd) )
        ? tau1 : tau2;

      // return p_k
      SLEQP_CALL(sleqp_sparse_vector_add_scaled(solver->z,
                                                solver->d,
                                                1.,
                                                tau,
                                                eps,
                                                newton_step));

      sleqp_log_debug("CG solver found negative curvature direction after %d iterations",
                      iteration);

      break;
    }

    // set alpha_j = (r_j^T * g_j) / (d_j^T * B_k * d_j)
    alpha = r_dot_g / dBd;

    // set z_{j+1} = z_j + alpha_j * d_j
    SLEQP_CALL(sleqp_sparse_vector_add_scaled(solver->z,
                                              solver->d,
                                              1.,
                                              alpha,
                                              eps,
                                              solver->sparse_cache));

    z_next_nrm_sq  = sleqp_sparse_vector_norm_sq(solver->sparse_cache);

    // if ||z_{j+1}|| >= Delta_k:

    if(z_next_nrm_sq >= trust_radius*trust_radius)
    {
      // find tau >= 0 such that p_k = z + tau*d_j satisfies ||p_k|| = Delta_k

      double z_dot_d;

      SLEQP_CALL(sleqp_sparse_vector_dot(solver->z, solver->d, &z_dot_d));

      const double d_nrm_sq = sleqp_sparse_vector_norm_sq(solver->d);
      const double z_nrm_sq = z_curr_nrm_sq;

      assert(d_nrm_sq > 0.);

      const double inner = z_dot_d*z_dot_d - d_nrm_sq*(z_nrm_sq - trust_radius*trust_radius);

      const double tau = 1./d_nrm_sq * (-z_dot_d + sqrt(inner));

      assert(tau >= 0);

      // set p_k = z_{j+1} + (tau-alpha)*d_j satisfies ||p_k|| = Delta_k
      // return p_k
      SLEQP_CALL(sleqp_sparse_vector_add_scaled(solver->z,
                                                solver->d,
                                                1.,
                                                tau,
                                                eps,
                                                newton_step));

      sleqp_num_assert(sleqp_is_eq(sleqp_sparse_vector_norm(newton_step),
                                   trust_radius,
                                   eps));

      sleqp_log_debug("CG solver found boundary solution after %d iterations",
                      iteration);

      break;
    }

    SLEQP_CALL(sleqp_sparse_vector_copy(solver->sparse_cache,
                                        solver->z));

    z_curr_nrm_sq = z_next_nrm_sq;

    // set r_{j+1} = r_j + alpha_j * (B_k * d_j)
    SLEQP_CALL(sleqp_sparse_vector_add_scaled(solver->r,
                                              solver->Bd,
                                              1.,
                                              alpha,
                                              eps,
                                              solver->sparse_cache));

    SLEQP_CALL(sleqp_sparse_vector_copy(solver->sparse_cache, solver->r));

    // set g_{j+1} = P[r_{j+1}]
    SLEQP_CALL(sleqp_aug_jacobian_projection(jacobian,
                                             solver->r,
                                             solver->g,
                                             NULL));

    // set beta_{j+1} = r_{j+1}^T * g_{j+1} / r_j^T * g_j
    beta = 1. / r_dot_g;
    SLEQP_CALL(sleqp_sparse_vector_dot(solver->r,
                                       solver->g,
                                       &r_dot_g));
    beta *= r_dot_g;

    // set d_{j+1} = - g_{j+1} + beta_{j+1} * d_j
    SLEQP_CALL(sleqp_sparse_vector_add_scaled(solver->g,
                                              solver->d,
                                              -1.,
                                              beta,
                                              eps,
                                              solver->sparse_cache));

    SLEQP_CALL(sleqp_sparse_vector_copy(solver->sparse_cache,
                                        solver->d));
  }
  // end loop

  sleqp_num_assert(sleqp_is_leq(sleqp_sparse_vector_norm(newton_step),
                                trust_radius,
                                eps));

  SLEQP_CALL(sleqp_timer_stop(solver->timer));

  return SLEQP_OKAY;
}


SLEQP_RETCODE
sleqp_steihaug_solver_create(SleqpTRSolver** solver_star,
                             SleqpProblem* problem,
                             SleqpParams* params,
                             SleqpOptions* options)
{
  SleqpSteihaugSolver *solver = NULL;

  SLEQP_CALL(sleqp_malloc(&solver));

  *solver = (SleqpSteihaugSolver) {0};

  solver->refcount = 1;

  solver->problem = problem;

  SLEQP_CALL(sleqp_params_capture(params));
  solver->params = params;

  solver->max_iter = sleqp_options_get_max_newton_iterations(options);

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->d, problem->num_variables));
  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->Bd, problem->num_variables));
  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->g, problem->num_variables));
  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->r, problem->num_variables));
  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->z, problem->num_variables));
  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->sparse_cache, problem->num_variables));

  SLEQP_CALL(sleqp_timer_create(&solver->timer));

  SleqpTRCallbacks callbacks = {
    .solve = steihaug_solver_solve,
    .free = steihaug_solver_free
  };

  SLEQP_CALL(sleqp_tr_solver_create(solver_star, problem, params, &callbacks, (void*) solver));

  return SLEQP_OKAY;
}
