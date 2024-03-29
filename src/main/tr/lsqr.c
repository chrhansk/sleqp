#include "lsqr.h"

#include <math.h>

#include "cmp.h"
#include "fail.h"
#include "mem.h"
#include "timer.h"
#include "tr_util.h"

struct SleqpLSQRSolver
{
  int refcount;

  double time_limit;

  SleqpSettings* settings;
  SleqpTimer* timer;

  int forward_dim;
  int adjoint_dim;

  SleqpLSQRCallbacks callbacks;
  void* data;

  SleqpVec* u;
  SleqpVec* v;
  SleqpVec* w;

  SleqpVec* p;
  SleqpVec* q;

  SleqpVec* d;

  SleqpVec* sparse_cache;
};

SLEQP_RETCODE
sleqp_lsqr_solver_create(SleqpLSQRSolver** star,
                         SleqpSettings* settings,
                         int forward_dim,
                         int adjoint_dim,
                         SleqpLSQRCallbacks* callbacks,
                         void* data)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpLSQRSolver* solver = *star;

  *solver = (SleqpLSQRSolver){0};

  solver->refcount   = 1;
  solver->time_limit = SLEQP_NONE;

  SLEQP_CALL(sleqp_settings_capture(settings));
  solver->settings = settings;

  SLEQP_CALL(sleqp_timer_create(&solver->timer));

  solver->forward_dim = forward_dim;
  solver->adjoint_dim = adjoint_dim;

  solver->callbacks = *callbacks;
  solver->data      = data;

  SLEQP_CALL(sleqp_vec_create_empty(&solver->u, 0));

  SLEQP_CALL(sleqp_vec_create_empty(&solver->v, solver->forward_dim));

  SLEQP_CALL(sleqp_vec_create_empty(&solver->w, solver->forward_dim));

  SLEQP_CALL(sleqp_vec_create_empty(&solver->p, 0));

  SLEQP_CALL(sleqp_vec_create_empty(&solver->q, solver->forward_dim));

  SLEQP_CALL(sleqp_vec_create_empty(&solver->d, solver->forward_dim));

  SLEQP_CALL(sleqp_vec_create_empty(&solver->sparse_cache, 0));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_lsqr_set_time_limit(SleqpLSQRSolver* solver, double time_limit)
{
  solver->time_limit = time_limit;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_lsqr_solver_resize(SleqpLSQRSolver* solver,
                         int forward_dim,
                         int adjoint_dim)
{
  solver->forward_dim = forward_dim;
  solver->adjoint_dim = adjoint_dim;

  SLEQP_CALL(sleqp_vec_resize(solver->v, forward_dim));
  SLEQP_CALL(sleqp_vec_resize(solver->w, forward_dim));

  SLEQP_CALL(sleqp_vec_resize(solver->q, forward_dim));
  SLEQP_CALL(sleqp_vec_resize(solver->d, forward_dim));

  SLEQP_CALL(sleqp_vec_resize(solver->u, adjoint_dim));
  SLEQP_CALL(sleqp_vec_resize(solver->p, adjoint_dim));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
normalize(SleqpVec* vec, double* norm)
{
  *norm = sleqp_vec_norm(vec);

  if (*norm)
  {
    SLEQP_CALL(sleqp_vec_scale(vec, 1. / (*norm)));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
forward_product(SleqpLSQRSolver* solver,
                const SleqpVec* direction,
                SleqpVec* product)
{
  assert(direction->dim == solver->forward_dim);
  assert(product->dim == solver->adjoint_dim);

  return solver->callbacks.prod_forward(direction, product, solver->data);
}

static SLEQP_RETCODE
adjoint_product(SleqpLSQRSolver* solver,
                const SleqpVec* direction,
                SleqpVec* product)
{
  assert(direction->dim == solver->adjoint_dim);
  assert(product->dim == solver->forward_dim);

  return solver->callbacks.prod_adjoint(direction, product, solver->data);
}

static SLEQP_RETCODE
compute_objective(SleqpLSQRSolver* solver,
                  const SleqpVec* rhs,
                  const SleqpVec* sol,
                  double* objective,
                  double* opt_res)
{
  const double zero_eps
    = sleqp_settings_real_value(solver->settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  SLEQP_CALL(forward_product(solver, sol, solver->p));

  SLEQP_CALL(
    sleqp_vec_add_scaled(solver->p, rhs, 1., -1., zero_eps, solver->u));

  SLEQP_CALL(adjoint_product(solver, solver->u, solver->w));

  const double res = sleqp_vec_norm(solver->u);

  (*objective) = .5 * (res * res);

  (*opt_res) = sleqp_vec_norm(solver->w);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_lsqr_solver_solve(SleqpLSQRSolver* solver,
                        const SleqpVec* rhs,
                        double rel_tol,
                        double trust_radius,
                        SleqpVec* sol)
{
  const int forward_dim = solver->forward_dim;
  int iteration         = 0;

  assert(rhs->dim == solver->adjoint_dim);
  assert(sol->dim == solver->forward_dim);

  const double eps = sleqp_settings_real_value(solver->settings, SLEQP_SETTINGS_REAL_EPS);

  const double zero_eps
    = sleqp_settings_real_value(solver->settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  sleqp_log_debug("Solving a least-squares subproblem with %d rows, %d columns",
                  solver->adjoint_dim,
                  solver->forward_dim);

  SLEQP_CALL(sleqp_timer_start(solver->timer));

  SleqpVec* x = sol;
  SleqpVec* u = solver->u;
  SleqpVec* v = solver->v;
  SleqpVec* w = solver->w;

  SleqpVec* p = solver->p;
  SleqpVec* q = solver->q;

  SleqpVec* d = solver->d;

  const SleqpVec* b = rhs;
  SleqpVec* t       = solver->sparse_cache;

  double alpha, beta;

  SLEQP_CALL(sleqp_vec_copy(b, u));
  SLEQP_CALL(normalize(u, &beta));

  SLEQP_CALL(adjoint_product(solver, u, v));
  SLEQP_CALL(normalize(v, &alpha));

  SLEQP_CALL(sleqp_vec_copy(v, w));

  SLEQP_CALL(sleqp_vec_clear(x));

  double phib = beta;
  double rhob = alpha;

  {
    const double res       = phib;
    const double objective = .5 * (res * res);
    const double opt_res   = phib * alpha;

    sleqp_log_debug("Initial objective: %e, residuum: %e", objective, opt_res);
  }

  bool reached_time_limit = false;

  for (iteration = 1; iteration <= forward_dim; ++iteration)
  {
    // Continue bidiagonalization
    SLEQP_CALL(forward_product(solver, v, p));
    SLEQP_CALL(sleqp_vec_add_scaled(p, u, 1., -alpha, zero_eps, t));
    SLEQP_CALL(sleqp_vec_copy(t, u));
    SLEQP_CALL(normalize(u, &beta));

    SLEQP_CALL(adjoint_product(solver, u, q));
    SLEQP_CALL(sleqp_vec_add_scaled(q, v, 1., -beta, zero_eps, t));
    SLEQP_CALL(sleqp_vec_copy(t, v));
    SLEQP_CALL(normalize(v, &alpha));

    // Construct next orthogonal transformation (Givens rotation)
    double rho     = hypot(rhob, beta);
    const double c = rhob / rho;
    const double s = beta / rho;
    double theta   = s * alpha;
    rhob           = (-c) * alpha;
    double phi     = c * phib;
    phib           = s * phib;

    SLEQP_CALL(sleqp_vec_add_scaled(x, w, 1., phi / rho, zero_eps, t));

    const double norm = sleqp_vec_norm(t);

    if (trust_radius != SLEQP_NONE && sleqp_is_gt(norm, trust_radius, eps))
    {
      SLEQP_CALL(sleqp_vec_add_scaled(x, t, -1., 1., zero_eps, d));

      SLEQP_CALL(
        sleqp_tr_compute_bdry_sol(x, d, solver->settings, trust_radius, t));

      SLEQP_CALL(sleqp_vec_copy(t, x));

#if SLEQP_DEBUG
      {
        double objective;
        double opt_res;

        SLEQP_CALL(compute_objective(solver, rhs, x, &objective, &opt_res));

        sleqp_log_debug("Iteration %d, objective: %e, residuum: %e",
                        iteration,
                        objective,
                        opt_res);
      }
#endif

      SLEQP_CALL(sleqp_timer_stop(solver->timer));

      sleqp_log_debug(
        "LSQR solver terminated with a boundary solution after %d iterations",
        iteration);

      return SLEQP_OKAY;
    }

    SLEQP_CALL(sleqp_vec_copy(t, x));

    SLEQP_CALL(sleqp_vec_add_scaled(v, w, 1., -theta / rho, zero_eps, t));

    SLEQP_CALL(sleqp_vec_copy(t, w));

    const double res       = phib;
    const double objective = .5 * (res * res);

    const double opt_res = phib * alpha * fabs(c);

    sleqp_log_debug("Iteration %d, objective: %e, residuum: %e",
                    iteration,
                    objective,
                    opt_res);

    if (opt_res <= rel_tol)
    {
      break;
    }

    if (solver->time_limit != SLEQP_NONE
        && sleqp_timer_elapsed(solver->timer) >= solver->time_limit)
    {
      reached_time_limit = true;
      break;
    }
  }

  SLEQP_CALL(sleqp_timer_stop(solver->timer));

  if (reached_time_limit)
  {
    return SLEQP_ABORT_TIME;
  }

  sleqp_log_debug(
    "LSQR solver terminated with an interior solution after %d iterations",
    iteration);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_lsqr_solver_capture(SleqpLSQRSolver* solver)
{
  ++solver->refcount;
  return SLEQP_OKAY;
}

static SLEQP_RETCODE
lsqr_solver_free(SleqpLSQRSolver** star)
{
  SleqpLSQRSolver* solver = *star;

  if (!solver)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_vec_free(&solver->sparse_cache));

  SLEQP_CALL(sleqp_vec_free(&solver->d));

  SLEQP_CALL(sleqp_vec_free(&solver->q));
  SLEQP_CALL(sleqp_vec_free(&solver->p));

  SLEQP_CALL(sleqp_vec_free(&solver->w));

  SLEQP_CALL(sleqp_vec_free(&solver->v));
  SLEQP_CALL(sleqp_vec_free(&solver->u));

  SLEQP_CALL(sleqp_timer_free(&solver->timer));

  SLEQP_CALL(sleqp_settings_release(&solver->settings));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_lsqr_solver_release(SleqpLSQRSolver** star)
{
  SleqpLSQRSolver* solver = *star;

  if (!solver)
  {
    return SLEQP_OKAY;
  }

  if (--(solver->refcount) == 0)
  {
    SLEQP_CALL(lsqr_solver_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
