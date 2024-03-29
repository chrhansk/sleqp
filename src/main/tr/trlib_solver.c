#include "tr_solver.h"

#include <math.h>
#include <string.h>

#include <trlib.h>

#include "cmp.h"
#include "error.h"
#include "fail.h"
#include "mem.h"

#include "sparse/mat.h"

static const double tolerance_factor = 1e-2;

typedef struct
{
  SleqpProblem* problem;
  SleqpSettings* settings;

  // trlib-related data:
  trlib_int_t trlib_maxiter;
  trlib_int_t trlib_h_pointer;
  trlib_int_t* trlib_iwork;
  trlib_flt_t* trlib_fwork;
  trlib_int_t* trlib_timinig;

  trlib_int_t iwork_size, fwork_size;

  SleqpVec* s;
  SleqpVec* g;
  SleqpVec* gm;

  SleqpVec* h;

  SleqpVec* v;
  SleqpVec* p;
  SleqpVec* Hp;

  SleqpVec* l;
  SleqpVec* h_lhs;
  SleqpVec* h_rhs;

  SleqpMat* Q;

  double* dense_cache;
  SleqpVec* sparse_cache;

  SleqpTimer* timer;
} SolverData;

static SLEQP_RETCODE
trlib_free(void** star)
{
  SolverData* data = (SolverData*)*star;

  if (!data)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_timer_free(&data->timer));

  SLEQP_CALL(sleqp_mat_release(&data->Q));

  SLEQP_CALL(sleqp_vec_free(&data->sparse_cache));
  sleqp_free(&data->dense_cache);

  SLEQP_CALL(sleqp_vec_free(&data->h_rhs));
  SLEQP_CALL(sleqp_vec_free(&data->h_lhs));

  SLEQP_CALL(sleqp_vec_free(&data->l));

  SLEQP_CALL(sleqp_vec_free(&data->Hp));
  SLEQP_CALL(sleqp_vec_free(&data->p));
  SLEQP_CALL(sleqp_vec_free(&data->v));

  SLEQP_CALL(sleqp_vec_free(&data->h));

  SLEQP_CALL(sleqp_vec_free(&data->gm));
  SLEQP_CALL(sleqp_vec_free(&data->g));
  SLEQP_CALL(sleqp_vec_free(&data->s));

  sleqp_free(&data->trlib_timinig);
  sleqp_free(&data->trlib_fwork);
  sleqp_free(&data->trlib_iwork);

  SLEQP_CALL(sleqp_settings_release(&data->settings));

  SLEQP_CALL(sleqp_problem_release(&data->problem));

  sleqp_free(star);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
matrix_push_column(SleqpMat* matrix, SleqpVec* vector, double scale)
{
  const int nnz = sleqp_mat_nnz(matrix);

  SLEQP_CALL(sleqp_mat_reserve(matrix, nnz + vector->nnz));

  const int num_cols = sleqp_mat_num_cols(matrix);
  const int num_rows = sleqp_mat_num_rows(matrix);

  const int column = num_cols;

  SLEQP_CALL(sleqp_mat_resize(matrix, num_rows, num_cols + 1));

  for (int k = 0; k < vector->nnz; ++k)
  {
    SLEQP_CALL(sleqp_mat_push(matrix,
                              vector->indices[k],
                              column,
                              vector->data[k] * scale));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
trlib_get_status_string(int value, const char** message)
{
  switch (value)
  {
  case TRLIB_CLR_CONTINUE:
    (*message) = "TRLIB_CLR_CONTINUE";
    break;
  case TRLIB_CLR_CONV_BOUND:
    (*message) = "TRLIB_CLR_CONV_BOUND";
    break;
  case TRLIB_CLR_CONV_INTERIOR:
    (*message) = "TRLIB_CLR_CONV_INTERIOR";
    break;
  case TRLIB_CLR_APPROX_HARD:
    (*message) = "TRLIB_CLR_APPROX_HARD";
    break;
  case TRLIB_CLR_NEWTON_BREAK:
    (*message) = "TRLIB_CLR_NEWTON_BREAK";
    break;
  case TRLIB_TTR_HARD_INIT_LAM:
    (*message) = "TRLIB_TTR_HARD_INIT_LAM";
    break;
  case TRLIB_CLR_ITMAX:
    (*message) = "TRLIB_CLR_ITMAX";
    break;
  case TRLIB_CLR_UNBDBEL:
    (*message) = "TRLIB_CLR_UNBDBEL";
    break;
  case TRLIB_CLR_FAIL_FACTOR:
    (*message) = "TRLIB_CLR_FAIL_FACTOR";
    break;
  case TRLIB_CLR_FAIL_LINSOLVE:
    (*message) = "TRLIB_CLR_FAIL_LINSOLVE";
    break;
  case TRLIB_CLR_FAIL_NUMERIC:
    (*message) = "TRLIB_CLR_FAIL_NUMERIC";
    break;
  case TRLIB_CLR_UNLIKE_CONV:
    (*message) = "TRLIB_CLR_UNLIKE_CONV";
    break;
  case TRLIB_CLR_PCINDEF:
    (*message) = "TRLIB_CLR_PCINDEF";
    break;
  case TRLIB_CLR_UNEXPECT_INT:
    (*message) = "TRLIB_CLR_UNEXPECT_INT";
    break;
  case TRLIB_CLR_FAIL_TTR:
    (*message) = "TRLIB_CLR_FAIL_TTR";
    break;
  case TRLIB_CLR_FAIL_HARD:
    (*message) = "TRLIB_CLR_FAIL_HARD";
    break;
  default:
    (*message) = "(unknown)";
    break;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
check_optimality(SolverData* data,
                 SleqpAugJac* jacobian,
                 const SleqpVec* multipliers,
                 const SleqpVec* gradient,
                 const SleqpVec* newton_step,
                 bool bdry_sol,
                 bool* is_optimal)
{
  SleqpProblem* problem = data->problem;

  const double eps = sleqp_settings_real_value(data->settings, SLEQP_SETTINGS_REAL_STAT_TOL);

  const double zero_eps
    = sleqp_settings_real_value(data->settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  (*is_optimal) = true;

  SLEQP_CALL(
    sleqp_problem_hess_prod(problem, newton_step, multipliers, data->Hp));

  SLEQP_CALL(sleqp_vec_add(gradient, data->Hp, zero_eps, data->g));

  SLEQP_CALL(sleqp_aug_jac_project_nullspace(jacobian, data->g, data->v));

  const double vnorm = sleqp_vec_norm(data->v);

  if (bdry_sol)
  {
    if (sleqp_is_zero(vnorm, eps))
    {
      return SLEQP_OKAY;
    }

    double dot;

    SLEQP_CALL(sleqp_vec_dot(data->v, newton_step, &dot));

    dot *= -1.;

    const double snorm = sleqp_vec_norm(newton_step);

    assert(snorm > 0.);

    if (!sleqp_is_eq(dot, vnorm * snorm, eps))
    {
      (*is_optimal) = false;
    }
  }
  else
  {
    if (!sleqp_is_zero(vnorm, eps))
    {
      (*is_optimal) = false;
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
trlib_loop(SolverData* data,
           SleqpAugJac* jacobian,
           const SleqpVec* multipliers,
           const SleqpVec* gradient,
           double trust_radius,
           double time_limit,
           trlib_int_t* trlib_ret)
{
  SleqpProblem* problem = data->problem;

  const int num_variables = sleqp_problem_num_vars(problem);

  const double inf = sleqp_infinity();

  const double zero_eps
    = sleqp_settings_real_value(data->settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  const double stat_eps
    = sleqp_settings_real_value(data->settings, SLEQP_SETTINGS_REAL_STAT_TOL);

  const double rel_tol = stat_eps * tolerance_factor;

  trlib_int_t equality      = 0;
  trlib_int_t maxlanczos    = data->trlib_maxiter;
  trlib_int_t ctl_invariant = 0;
  trlib_int_t refine        = 1;
  trlib_int_t verbose       = (sleqp_log_level() >= SLEQP_LOG_DEBUG);
  trlib_int_t unicode       = 0;
  trlib_flt_t tol_rel_i     = rel_tol;
  trlib_flt_t tol_abs_i     = zero_eps;
  trlib_flt_t tol_rel_b     = rel_tol;
  trlib_flt_t tol_abs_b     = zero_eps;
  trlib_flt_t obj_lo        = -inf;
  trlib_int_t convexify     = 1;
  trlib_int_t earlyterm     = 1;

  trlib_int_t init = 0;

  trlib_flt_t v_dot_g = 0.0, p_dot_Hp = 0.0, g_dot_g = 0.0, flt1, flt2, flt3;

  trlib_int_t action, ityp;

  // trlib_flt_t zero = TRLIB_EPS*TRLIB_EPS;

  trlib_flt_t zero = zero_eps;

  trlib_int_t iter = 0;

  trlib_int_t* iwork = data->trlib_iwork;
  trlib_flt_t* fwork = data->trlib_fwork;

  trlib_int_t* timing = data->trlib_timinig;

  trlib_int_t maxiter = data->trlib_maxiter;

  SLEQP_CALL(sleqp_timer_start(data->timer));

  {
    init = TRLIB_CLS_INIT;
    memset(data->trlib_fwork, 0, (data->fwork_size) * sizeof(trlib_flt_t));
    memset(data->trlib_iwork, 0, (data->iwork_size) * sizeof(trlib_int_t));

    trlib_krylov_prepare_memory(maxiter, fwork);
  }

  SLEQP_CALL(sleqp_mat_clear(data->Q));

  SLEQP_CALL(sleqp_vec_clear(data->p));
  SLEQP_CALL(sleqp_vec_clear(data->Hp));

  char* prefix = "trlib: ";
  FILE* fout   = stderr;

  bool exhausted_time_limit = false;

  while (1)
  {
    (*trlib_ret) = trlib_krylov_min(init,          // trlib_int_t init
                                    trust_radius,  // trlib_flt_t radius
                                    equality,      // trlib_int_t equality
                                    maxiter,       // trlib_int_t itmax
                                    maxlanczos,    // trlib_int_t itmax_lanczos
                                    tol_rel_i,     // trlib_flt_t tol_rel_i
                                    tol_abs_i,     // trlib_flt_t tol_abs_i
                                    tol_rel_b,     // trlib_flt_t tol_rel_b
                                    tol_abs_b,     // trlib_flt_t tol_abs_b
                                    zero,          // trlib_flt_t zero
                                    obj_lo,        // trlib_flt_t obj_lo
                                    ctl_invariant, // trlib_int_t ctl_invariant
                                    convexify,     // trlib_int_t convexify
                                    earlyterm,     // trlib_int_t earlyterm
                                    g_dot_g,       // trlib_flt_t g_dot_g
                                    v_dot_g,       // trlib_flt_t v_dot_g
                                    p_dot_Hp,      // trlib_flt_t p_dot_Hp
                                    iwork,         // trlib_int_t *iwork
                                    fwork,         // trlib_flt_t *fwork
                                    refine,        // trlib_int_t refine
                                    verbose,       // trlib_int_t verbose
                                    unicode,       // trlib_int_t unicode
                                    prefix,        // char *prefix
                                    fout,          // FILE *fout
                                    timing,        // trlib_int_t *timing
                                    &action,       // trlib_int_t *action
                                    &iter,         // trlib_int_t *iter
                                    &ityp,         // trlib_int_t *ityp
                                    &flt1,         // trlib_flt_t *flt1
                                    &flt2,         // trlib_flt_t *flt2
                                    &flt3);        // trlib_flt_t *flt3

    init = 0;

    switch (action)
    {
    case TRLIB_CLA_TRIVIAL:
    {
      break;
    }
    case TRLIB_CLA_INIT:
    {
      // reset s to 0
      SLEQP_CALL(sleqp_vec_clear(data->s));

      SLEQP_CALL(sleqp_vec_copy(gradient, data->g));

      SLEQP_CALL(sleqp_vec_clear(data->gm));

      SLEQP_CALL(sleqp_aug_jac_project_nullspace(jacobian, data->g, data->v));

      SLEQP_CALL(sleqp_vec_copy(data->v, data->p));

      SLEQP_CALL(sleqp_vec_scale(data->p, -1.));

      g_dot_g = sleqp_vec_norm_sq(data->g);

      SLEQP_CALL(sleqp_vec_dot(data->v, data->g, &v_dot_g));

      SLEQP_CALL(
        sleqp_problem_hess_prod(problem, data->p, multipliers, data->Hp));

      SLEQP_CALL(sleqp_vec_dot(data->p, data->Hp, &p_dot_Hp));

      const int num_rows = sleqp_mat_num_rows(data->Q);

      SLEQP_CALL(sleqp_mat_resize(data->Q, num_rows, 0));
      SLEQP_CALL(sleqp_mat_clear(data->Q));

      assert(sleqp_mat_num_cols(data->Q) == 0);

      // assert(v_dot_g > 0);

      double scale = sqrt(v_dot_g);

      if (sleqp_is_pos(scale, zero_eps))
      {
        scale = 1. / scale;

        SLEQP_CALL(matrix_push_column(data->Q, data->v, scale));

        assert(sleqp_mat_is_valid(data->Q));
      }

      break;
    }
    case TRLIB_CLA_RETRANSF:
    {
      SLEQP_CALL(sleqp_vec_set_from_raw(data->h,
                                        fwork + data->trlib_h_pointer,
                                        iter + 1,
                                        zero_eps));

      SLEQP_CALL(sleqp_vec_resize(data->h, sleqp_mat_num_cols(data->Q)));

      SLEQP_CALL(sleqp_mat_mult_vec(data->Q, data->h, data->dense_cache));

      SLEQP_CALL(sleqp_vec_set_from_raw(data->s,
                                        data->dense_cache,
                                        num_variables,
                                        zero_eps));

      break;
    }
    case TRLIB_CLA_UPDATE_STATIO:
    {
      if (ityp == TRLIB_CLT_CG)
      {
        SLEQP_CALL(sleqp_vec_add_scaled(data->s,
                                        data->p,
                                        1.,
                                        flt1,
                                        zero_eps,
                                        data->sparse_cache));

        SLEQP_CALL(sleqp_vec_copy(data->sparse_cache, data->s));
      }
      break;
    }
    case TRLIB_CLA_UPDATE_GRAD:
    {
      if (ityp == TRLIB_CLT_CG)
      {
        if (iter + 1 == sleqp_mat_num_cols(data->Q))
        {
          const int num_rows = sleqp_mat_num_rows(data->Q);
          SLEQP_CALL(sleqp_mat_resize(data->Q, num_rows, iter));
        }

        assert(iter == sleqp_mat_num_cols(data->Q));

        SLEQP_CALL(matrix_push_column(data->Q, data->v, flt2));

        assert(iter + 1 == sleqp_mat_num_cols(data->Q));
        assert(sleqp_mat_is_valid(data->Q));

        SLEQP_CALL(sleqp_vec_copy(data->g, data->gm));

        SLEQP_CALL(sleqp_vec_add_scaled(data->g,
                                        data->Hp,
                                        1.,
                                        flt1,
                                        zero_eps,
                                        data->sparse_cache));

        SLEQP_CALL(sleqp_vec_copy(data->sparse_cache, data->g));

        SLEQP_CALL(sleqp_aug_jac_project_nullspace(jacobian, data->g, data->v));

        g_dot_g = sleqp_vec_norm_sq(data->g);

        SLEQP_CALL(sleqp_vec_dot(data->v, data->g, &v_dot_g));
      }
      if (ityp == TRLIB_CLT_L)
      {
        SLEQP_CALL(sleqp_vec_add_scaled(data->Hp,
                                        data->g,
                                        1.,
                                        flt1,
                                        zero_eps,
                                        data->sparse_cache));

        SLEQP_CALL(sleqp_vec_add_scaled(data->sparse_cache,
                                        data->gm,
                                        1.,
                                        flt2,
                                        zero_eps,
                                        data->s));

        SLEQP_CALL(sleqp_vec_copy(data->g, data->gm));

        SLEQP_CALL(sleqp_vec_scale(data->gm, flt3));

        SLEQP_CALL(sleqp_vec_copy(data->s, data->g));

        SLEQP_CALL(sleqp_aug_jac_project_nullspace(jacobian, data->g, data->v));

        SLEQP_CALL(sleqp_vec_dot(data->v, data->g, &v_dot_g));
      }

      break;
    }
    case TRLIB_CLA_CONV_HARD:
    {
      SLEQP_CALL(sleqp_problem_hess_prod(problem,
                                         data->s,
                                         multipliers,
                                         data->sparse_cache));

      SLEQP_CALL(
        sleqp_vec_add(data->sparse_cache, gradient, zero_eps, data->l));

      SLEQP_CALL(sleqp_vec_add_scaled(data->l,
                                      data->s,
                                      1.,
                                      flt1,
                                      zero_eps,
                                      data->h_lhs));

      SLEQP_CALL(
        sleqp_aug_jac_project_nullspace(jacobian, data->l, data->sparse_cache));

      SLEQP_CALL(sleqp_vec_add_scaled(data->sparse_cache,
                                      data->s,
                                      1.,
                                      flt1,
                                      zero_eps,
                                      data->h_rhs));

      SLEQP_CALL(sleqp_vec_dot(data->h_lhs, data->h_rhs, &v_dot_g));

      break;
    }

    case TRLIB_CLA_NEW_KRYLOV:
    {
      // This should not happen in our case...
      assert(false);
      break;
    }
    case TRLIB_CLA_UPDATE_DIR:
    {
      if (ityp == TRLIB_CLT_CG)
      {
        assert(flt1 == -1.);

        SLEQP_CALL(sleqp_vec_add_scaled(data->v,
                                        data->p,
                                        flt1,
                                        flt2,
                                        zero_eps,
                                        data->sparse_cache));

        SLEQP_CALL(sleqp_vec_copy(data->sparse_cache, data->p));

        SLEQP_CALL(sleqp_problem_hess_prod(problem,
                                           data->p,
                                           multipliers,
                                           data->Hp));

        SLEQP_CALL(sleqp_vec_dot(data->p, data->Hp, &p_dot_Hp));
      }
      else if (ityp == TRLIB_CLT_L)
      {
        assert(flt2 == 0.);

        const int num_cols = sleqp_mat_num_cols(data->Q);

        if (iter + 1 == num_cols)
        {
          const int num_rows = sleqp_mat_num_rows(data->Q);
          SLEQP_CALL(sleqp_mat_resize(data->Q, num_rows, iter));
          assert(sleqp_mat_is_valid(data->Q));
        }

        assert(iter == sleqp_mat_num_cols(data->Q));

        SLEQP_CALL(sleqp_vec_add_scaled(data->v,
                                        data->p,
                                        flt1,
                                        flt2,
                                        zero_eps,
                                        data->sparse_cache));

        SLEQP_CALL(sleqp_vec_copy(data->sparse_cache, data->p));

        SLEQP_CALL(sleqp_problem_hess_prod(problem,
                                           data->p,
                                           multipliers,
                                           data->Hp));

        SLEQP_CALL(sleqp_vec_dot(data->p, data->Hp, &p_dot_Hp));

        SLEQP_CALL(matrix_push_column(data->Q, data->p, 1.));

        assert(iter + 1 == sleqp_mat_num_cols(data->Q));

        assert(sleqp_mat_is_valid(data->Q));
      }
      break;
    }
    case TRLIB_CLA_OBJVAL:
    {
      double lin_term, quad_term;

      SLEQP_CALL(sleqp_problem_hess_bilinear(problem,
                                             data->s,
                                             multipliers,
                                             &quad_term));

      SLEQP_CALL(sleqp_vec_dot(data->s, data->g, &lin_term));

      g_dot_g = lin_term + .5 * quad_term;

      break;
    }
    default:
    {
      sleqp_raise(SLEQP_INTERNAL_ERROR,
                  "Invalid trlib action requested: %ld",
                  action);
    }
    }

    if ((*trlib_ret) < TRLIB_CLR_CONTINUE)
    {
      break;
    }

    if (time_limit != SLEQP_NONE
        && sleqp_timer_elapsed(data->timer) >= time_limit)
    {
      exhausted_time_limit = true;
      break;
    }
  }

  SLEQP_CALL(sleqp_timer_stop(data->timer));

  if (exhausted_time_limit)
  {
    return SLEQP_ABORT_TIME;
  }
  else
  {
    assert((*trlib_ret) != TRLIB_CLR_CONTINUE);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
trlib_rayleigh(double* min_rayleigh, double* max_rayleigh, void* solver_data)
{
  SolverData* data = (SolverData*)solver_data;

  (*min_rayleigh) = data->trlib_fwork[14];
  (*max_rayleigh) = data->trlib_fwork[13];

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
trlib_solve(SleqpAugJac* jacobian,
            const SleqpVec* multipliers,
            const SleqpVec* gradient,
            SleqpVec* newton_step,
            double trust_radius,
            double* tr_dual,
            double time_limit,
            void* solver_data)
{
  trlib_int_t ret = 0;

  SolverData* data = (SolverData*)solver_data;

  SLEQP_CALL(sleqp_vec_clear(newton_step));

  SLEQP_CALL(trlib_loop(data,
                        jacobian,
                        multipliers,
                        gradient,
                        trust_radius,
                        time_limit,
                        &ret));

  // trust region dual, located int trlib
  *tr_dual = data->trlib_fwork[7];

  // We may loose some orthogonality in the process, reproject to
  // be sure
  SLEQP_CALL(sleqp_aug_jac_project_nullspace(jacobian, data->s, newton_step));

  // check optimality

  const bool solved_optimally
    = (ret == TRLIB_CLR_CONV_BOUND) || (ret == TRLIB_CLR_CONV_INTERIOR);

  const bool converged_bdry = (ret == TRLIB_CLR_CONV_BOUND);

#if SLEQP_DEBUG
  if (solved_optimally)
  {
    bool tr_subproblem_is_optimal;

    SLEQP_CALL(check_optimality(data,
                                jacobian,
                                multipliers,
                                gradient,
                                newton_step,
                                converged_bdry,
                                &tr_subproblem_is_optimal));

    sleqp_num_assert(tr_subproblem_is_optimal);
  }

  // TODO: Choose appropriate tolerance
  // assert(sleqp_is_leq(sleqp_vec_norm(newton_step), trust_radius,
  // eps));

#endif

  const char* trlib_status_string;

  SLEQP_CALL(trlib_get_status_string(ret, &trlib_status_string));

  if (!solved_optimally)
  {
    sleqp_log_warn("Failed to solve trust region subproblem, reason: %s",
                   trlib_status_string);
  }

  const bool failure_occured
    = (ret == TRLIB_CLR_FAIL_FACTOR) || (ret == TRLIB_CLR_FAIL_LINSOLVE)
      || (ret == TRLIB_CLR_FAIL_NUMERIC) || (ret == TRLIB_CLR_FAIL_TTR);

  if (failure_occured)
  {
    sleqp_raise(SLEQP_INTERNAL_ERROR,
                "Failure %s occured in trlib",
                trlib_status_string);
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_trlib_solver_create(SleqpTRSolver** solver_star,
                          SleqpProblem* problem,
                          SleqpSettings* settings)
{
  SolverData* data = NULL;

  const int num_constraints = sleqp_problem_num_cons(problem);
  const int num_variables   = sleqp_problem_num_vars(problem);

  SLEQP_CALL(sleqp_malloc(&data));

  *data = (SolverData){0};

  data->problem = problem;
  SLEQP_CALL(sleqp_problem_capture(data->problem));

  SLEQP_CALL(sleqp_settings_capture(settings));
  data->settings = settings;

  const int max_newton_iter
    = sleqp_settings_int_value(settings, SLEQP_SETTINGS_INT_MAX_NEWTON_ITERATIONS);

  data->trlib_maxiter = num_variables;

  if (max_newton_iter != SLEQP_NONE)
  {
    data->trlib_maxiter = SLEQP_MIN(data->trlib_maxiter, max_newton_iter);
  }

  trlib_krylov_memory_size(data->trlib_maxiter,
                           &data->iwork_size,
                           &data->fwork_size,
                           &data->trlib_h_pointer);

  SLEQP_CALL(sleqp_alloc_array(&data->trlib_iwork, data->iwork_size));
  SLEQP_CALL(sleqp_alloc_array(&data->trlib_fwork, data->fwork_size));
  SLEQP_CALL(
    sleqp_alloc_array(&data->trlib_timinig, trlib_krylov_timing_size()));

  SLEQP_CALL(sleqp_vec_create_empty(&data->s, num_variables));
  SLEQP_CALL(sleqp_vec_create_empty(&data->g, num_variables));
  SLEQP_CALL(sleqp_vec_create_empty(&data->gm, num_variables));

  SLEQP_CALL(sleqp_vec_create_empty(&data->h, num_variables));

  SLEQP_CALL(sleqp_vec_create_empty(&data->v, num_variables));
  SLEQP_CALL(sleqp_vec_create_empty(&data->p, num_variables));
  SLEQP_CALL(sleqp_vec_create_empty(&data->Hp, num_variables));

  SLEQP_CALL(sleqp_vec_create_empty(&data->l, num_variables));

  SLEQP_CALL(sleqp_vec_create_empty(&data->h_lhs, num_variables));
  SLEQP_CALL(sleqp_vec_create_empty(&data->h_rhs, num_variables));

  SLEQP_CALL(
    sleqp_mat_create(&data->Q, num_variables, data->trlib_maxiter + 1, 0));

  SLEQP_CALL(sleqp_alloc_array(&data->dense_cache,
                               SLEQP_MAX(num_variables, num_constraints)));

  SLEQP_CALL(sleqp_vec_create_empty(&data->sparse_cache, num_variables));

  SLEQP_CALL(sleqp_timer_create(&data->timer));

  SleqpTRCallbacks callbacks
    = {.solve = trlib_solve, .rayleigh = trlib_rayleigh, .free = trlib_free};

  SLEQP_CALL(sleqp_tr_solver_create(solver_star, &callbacks, (void*)data));

  return SLEQP_OKAY;
}
