/*
 * An iterative solver for the trust-region subproblem using projected Conjugate Gradients (CG)
 * with Steihaug's modification for the boundary case. The augmented Jacobian system is used to 
 * project onto the nullspace of the active set identified in the LP step. The (1,1)-block of the 
 * projector is currently the identity, but could contain a Hessian preconditioner.
 */

#include "sleqp_tr_solver.h"

#include <math.h>
#include <trlib.h>

#include "sleqp_assert.h"
#include "sleqp_cmp.h"
#include "sleqp_mem.h"

typedef struct SleqpTRSolver SleqpTRSolver;


struct SleqpTRSolver
{
  int refcount;
  SleqpProblem* problem;
  SleqpParams* params;

  SleqpTimer* timer;
  double time_limit;

  SleqpSparseVec* d;            // projected descent direction
  SleqpSparseVec* Bd;           // Hessian-times-direction product
  SleqpSparseVec* g;            // projected residual
  SleqpSparseVec* r;            // full-space residual
  SleqpSparseVec* z;            // CG iterate

  SleqpSparseVec* sparse_cache; // a temporary t for computing t <- a*x + b*y; x <- t
};


SLEQP_RETCODE 
sleqp_tr_solver_create (
    SleqpTRSolver** star,
    SleqpProblem* problem,
    SleqpParams* params,
    SleqpOptions* options
)
{
    SLEQP_CALL(sleqp_malloc(star));

    SleqpTRSolver *solver = *star;

    *solver = (SleqpTRSolver) {0};

    solver->refcount = 1;

    solver->problem = problem;

    SLEQP_CALL (sleqp_params_capture (params));
    solver->params = params;

    solver->trlib_maxiter = problem->num_variables;
    const int max_newton_iter = sleqp_options_get_max_newton_iterations (options);
    if (max_newton_iter != SLEQP_NONE)
    {
        solver->trlib_maxiter = SLEQP_MIN (solver->trlib_maxiter, max_newton_iter);
    }

    SLEQP_CALL( sleqp_sparse_vector_create_empty (&solver->d, problem->num_variables));
    SLEQP_CALL( sleqp_sparse_vector_create_empty (&solver->Bd, problem->num_variables));
    SLEQP_CALL( sleqp_sparse_vector_create_empty (&solver->g, problem->num_variables));
    SLEQP_CALL( sleqp_sparse_vector_create_empty (&solver->r, problem->num_variables));
    SLEQP_CALL( sleqp_sparse_vector_create_empty (&solver->z, problem->num_variables));
    SLEQP_CALL( sleqp_sparse_vector_create_empty (&solver->sparse_cache, problem->num_variables));

    solver->time_limit = SLEQP_NONE;
    SLEQP_CALL (sleqp_timer_create (&data->timer));

    return SLEQP_OKAY;
}


static SLEQP_RETCODE 
tr_solver_free (
    SleqpTRSolver **star
)
{
    SleqpTRSolver* solver = *star;

    if (!solver)
    {
        return SLEQP_OKAY;
    }

    SLEQP_CALL (sleqp_timer_free (&solver->timer));

    SLEQP_CALL (sleqp_sparse_vector_free (&solver->sparse_cache));
    SLEQP_CALL (sleqp_sparse_vector_free (&solver->z));
    SLEQP_CALL (sleqp_sparse_vector_free (&solver->r));
    SLEQP_CALL (sleqp_sparse_vector_free (&solver->g));
    SLEQP_CALL (sleqp_sparse_vector_free (&solver->Bd));
    SLEQP_CALL (sleqp_sparse_vector_free (&solver->d));

    SLEQP_CALL (sleqp_params_release (&solver->params));

    sleqp_free (star);

    return SLEQP_OKAY;
}


SLEQP_RETCODE 
sleqp_tr_solver_set_time_limit (
    SleqpTRSolver* solver,
    double time_limit

{
    data->time_limit = time_limit;

    return SLEQP_OKAY;
}


SleqpTimer * 
sleqp_tr_solver_get_solve_timer (
    SleqpTRSolver* solver
)
{
    return data->timer;
}


SLEQP_RETCODE 
sleqp_tr_solver_capture (
    SleqpTRSolver* solver
)
{
    ++data->refcount;

    return SLEQP_OKAY;
}


SLEQP_RETCODE 
sleqp_tr_solver_release (
    SleqpTRSolver** star
)
{
    SleqpTRSolver* data = *star;

    if(!data)
    {
        return SLEQP_OKAY;
    }

    if(--(data->refcount) == 0)
    {
        SLEQP_CALL(tr_solver_free(star));
    }

    *star = NULL;

    return SLEQP_OKAY;
}


SLEQP_RETCODE 
sleqp_tr_solver_solve (
    SleqpTRSolver* solver,
    SleqpAugJacobian* jacobian,
    SleqpSparseVec* multipliers,
    SleqpSparseVec* gradient,
    SleqpSparseVec* newton_step,
    double trust_radius
)
{
    const long n = gradient->dim;
    const double one = 1.0;
    const double rel_tol = sleqp_params_get (data->params, SLEQP_PARAM_NEWTON_RELATIVE_TOL);
    double dBd;
    double alpha;
    double beta;
    double nrm_sq;
    double d_nrm_sq;
    double tau;
    double tau1;
    double tau2;
    double gd;
    double zBd;

    SLEQP_CALL (sleqp_timer_start (solver->timer));

    SLEQP_CALL(sleqp_sparse_vector_clear(newton_step));

    // set z0 such that P[z0] = 0
    sleqp_sparse_vector_clear (solver->z);

    // set r0 = nabla f_k
    sleqp_sparse_vector_copy (gradient, solver->r);

    // set g0 = P[r0]
    SLEQP_CALL (sleqp_aug_jacobian_projection (jacobian, solver->r, solver->g, NULL));

    // set d0 = -P[r0] = - P[nabla f_k] = -g0
    sleqp_sparse_vector_copy (solver->g, solver->d);
    sleqp_sparse_vector_scale (solver->d, -1.0);

    // if ||d0|| < eps_k: return p_k = P[z_0] = 0
    d_nrm_sq = sleqp_sparse_vector_norm_sq (solver->d);
    if (d_nrm_sq < rel_tol*rel_tol)
    {
        sleqp_sparse_vector_copy (solver->z, newton_step);
        SLEQP_CALL (sleqp_timer_stop (solver->timer));
        return SLEQP_OKAY;
    }

    // compute r_j^T * g_j
    r_dot_g = sleqp_sparse_vector_dot (solver->r, solver->g);

    // loop over pCG iterations j
    for (long j = 0; j < data->trlib_maxiter; ++j)
    {
        // if |r_{j+1}^T * g_{j+1}| < eps_k:
        if (abs (r_dot_g) < rel_tol)
        {
            // return p_k = z_{j+1}
            sleqp_sparse_vector_copy (solver->z, newton_step);
            break;
        }

        // compute B_k * d_j
        SLEQP_CALL (sleqp_func_hess_prod (func, &one, solver->d, multipliers, solver->Bd));

        // compute d_j^T * (B_k * d_j)
        sleqp_sparse_vector_dot (solver->d, solver->Bd, &dBd);

        // if d_j^T * B_k * d_j <= 0:
        if (dBd <= 0.0)
        {
            // z is a feasible iterate and d is a null-space direction, so P^T[p_k] = 0
            // find tau such that p_k = z_j + tau * d_j solves
            // min_{tau}   m_k(p) = f_k + nabla f_k^T * p + 0.5 * p^T * B_k * p 
            // subject to  ||p_k|| = Delta

            // solving for tau results in tau = - z^T * d / ||d||^2 \pm Delta / ||d|| 
            tau = sleqp_sparse_vector_dot (solver->z, solver->d);
            nrm_sq = sleqp_sparse_vector_dot (solver->d, solver->d);
            tau1 = - tau / nrm_sq + trust_radius / sqrt (nrm_sq);        /// @todo verify again
            tau2 = - tau / nrm_sq - trust_radius / sqrt (nrm_sq);   

            // compute dot products g^T * d and z^T * (B * d)
            gd = sleqp_sparse_vector_dot (gradient, solver->d);
            zBd = sleqp_sparse_vector_dot (solver->z, solver->Bd);

            // pick the boundary intersection with smaller quadratic model value
            // q(p) = tau * (g^T * d + z^T * B * d) + 0.5 * tau^2 * d^T * B * d + const.
            tau = (tau1*((gd + zBd) + 0.5*tau1*dBd) < tau2*((gd + zBd) + 0.5*tau2*dBd) ) 
                ? tau1 : tau2;

            // return p_k
            sleqp_sparse_vector_add_scaled (solver->z, solver->d, 1.0, tau, newton_step);
            break;
        }

        // set alpha_j = (r_j^T * g_j) / (d_j^T * B_k * d_j)
        alpha = r_dot_g / dBd;
        
        // set z_{j+1} = z_j + alpha_j * d_j
        sleqp_sparse_vector_add_scaled (solver->z, solver->d, 1.0, alpha, solver->sparse_cache);
        sleqp_sparse_vector_copy (solver->sparse_cache, solver->z);

        // if ||z_{j+1}|| >= Delta_k:
        nrm_sq = sleqp_sparse_vector_nrm_sq (solver->z);
        if (nrm_sq >= trust_radius*trust_radius)
        {
            // find tau >= 0 such that p_k = z + tau*d_j satisfies ||p_k|| = Delta_k
                
            // solving for tau results in tau = -z * d / ||d||^2 +/- Delta / ||d|| 
            tau = sleqp_sparse_vector_dot (solver->z, solver->d);
            nrm_sq = sleqp_sparse_vector_dot (solver->d, solver->d);
            tau = - tau / nrm_sq + trust_radius / sqrt (nrm_sq);        /// @todo verify again

            // set p_k = z_{j+1} + (tau-alpha)*d_j satisfies ||p_k|| = Delta_k
            // return p_k
            sleqp_sparse_vector_add_scaled (solver->z, solver->d, 1.0, tau-alpha, newton_step);
            break;
        }

        // set r_{j+1} = r_j + alpha_j * (B_k * d_j)
        sleqp_sparse_vector_add_scaled (solver->r, solver->Bd, 1.0, alpha, solver->sparse_cache);
        sleqp_sparse_vector_copy (solver->sparse_cache, solver->r);

        // set g_{j+1} = P[r_{j+1}]
        SLEQP_CALL (sleqp_aug_jacobian_projection (jacobian, solver->r, solver->g, NULL));

        // set beta_{j+1} = r_{j+1}^T * g_{j+1} / r_j^T * g_j
        beta = 1.0 / r_dot_g;
        r_dot_g = sleqp_sparse_vector_dot (solver->r, solver->g);
        beta *= r_dot_g;

        // set d_{j+1} = - g_{j+1} + beta_{j+1} * d_j
        sleqp_sparse_vector_add_scaled (solver->g, solver->d, -1.0, beta, solver->d);
    }
    // end loop

    SLEQP_CALL (sleqp_timer_stop (solver->timer));
    return SLEQP_OKAY;
}
