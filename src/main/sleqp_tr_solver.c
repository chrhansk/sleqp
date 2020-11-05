#include "sleqp_tr_solver.h"

#include <math.h>
#include <trlib.h>

#include "sleqp_cmp.h"
#include "sleqp_mem.h"

static const bool collect_rayleigh = false;

struct SleqpTRSolver
{
  int refcount;

  SleqpProblem* problem;
  SleqpParams* params;

  // trlib-related data:
  trlib_int_t trlib_maxiter;
  trlib_int_t trlib_h_pointer;
  trlib_int_t *trlib_iwork;
  trlib_flt_t *trlib_fwork;
  trlib_int_t *trlib_timinig;

  SleqpSparseVec* s;
  SleqpSparseVec* g;
  SleqpSparseVec* gm;

  SleqpSparseVec* h;

  SleqpSparseVec* v;
  SleqpSparseVec* p;
  SleqpSparseVec* Hp;

  SleqpSparseVec* l;
  SleqpSparseVec* h_lhs;
  SleqpSparseVec* h_rhs;

  SleqpSparseMatrix* Q;

  double* dense_cache;
  SleqpSparseVec* sparse_cache;

  SleqpTimer* timer;
  double time_limit;
};


SLEQP_RETCODE sleqp_tr_solver_create(SleqpTRSolver** star,
                                     SleqpProblem* problem,
                                     SleqpParams* params,
                                     SleqpOptions* options)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpTRSolver* data = *star;

  *data = (SleqpTRSolver) {0};

  data->refcount = 1;

  data->problem = problem;

  SLEQP_CALL(sleqp_params_capture(params));
  data->params = params;

  const int max_newton_iter = sleqp_options_get_max_newton_iterations(options);

  data->trlib_maxiter = problem->num_variables;

  if(max_newton_iter != SLEQP_NONE)
  {
    data->trlib_maxiter = SLEQP_MIN(data->trlib_maxiter, max_newton_iter);
  }

  trlib_int_t iwork_size, fwork_size;

  trlib_krylov_memory_size(data->trlib_maxiter,
                           &iwork_size,
                           &fwork_size,
                           &data->trlib_h_pointer);

  SLEQP_CALL(sleqp_calloc(&data->trlib_iwork, iwork_size));
  SLEQP_CALL(sleqp_calloc(&data->trlib_fwork, fwork_size));
  SLEQP_CALL(sleqp_calloc(&data->trlib_timinig, trlib_krylov_timing_size()));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->s,
                                              problem->num_variables));
  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->g,
                                              problem->num_variables));
  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->gm,
                                              problem->num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->h,
                                              problem->num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->v,
                                              problem->num_variables));
  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->p,
                                              problem->num_variables));
  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->Hp,
                                              problem->num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->l,
                                              problem->num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->h_lhs,
                                              problem->num_variables));
  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->h_rhs,
                                              problem->num_variables));

  SLEQP_CALL(sleqp_sparse_matrix_create(&data->Q,
                                        problem->num_variables,
                                        data->trlib_maxiter + 1,
                                        0));

  SLEQP_CALL(sleqp_calloc(&data->dense_cache,
                          SLEQP_MAX(problem->num_variables, problem->num_constraints)));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->sparse_cache,
                                              problem->num_variables));

  data->time_limit = SLEQP_NONE;

  SLEQP_CALL(sleqp_timer_create(&data->timer));
}

static SLEQP_RETCODE tr_solver_free(SleqpTRSolver** star)
{
  SleqpTRSolver* data = *star;

  if(!data)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_timer_free(&data->timer));

  SLEQP_CALL(sleqp_sparse_matrix_release(&data->Q));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->sparse_cache));
  sleqp_free(&data->dense_cache);

  SLEQP_CALL(sleqp_sparse_vector_free(&data->h_rhs));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->h_lhs));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->l));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->Hp));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->p));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->v));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->h));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->gm));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->g));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->s));

  sleqp_free(&data->trlib_timinig);
  sleqp_free(&data->trlib_fwork);
  sleqp_free(&data->trlib_iwork);

  SLEQP_CALL(sleqp_params_release(&data->params));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_tr_solver_set_time_limit(SleqpTRSolver* data,
                                             double time_limit)
{
  data->time_limit = time_limit;

  return SLEQP_OKAY;
}

SleqpTimer* sleqp_tr_solver_get_solve_timer(SleqpTRSolver* data)
{
  return data->timer;
}

static SLEQP_RETCODE matrix_push_column(SleqpSparseMatrix* matrix,
                                        SleqpSparseVec* vector,
                                        double scale)
{
  const int nnz = sleqp_sparse_matrix_get_nnz(matrix);

  SLEQP_CALL(sleqp_sparse_matrix_reserve(matrix, nnz + vector->nnz));

  const int num_cols = sleqp_sparse_matrix_get_num_cols(matrix);
  const int num_rows = sleqp_sparse_matrix_get_num_rows(matrix);

  const int column = num_cols;

  SLEQP_CALL(sleqp_sparse_matrix_resize(matrix,
                                        num_rows,
                                        num_cols + 1));

  for(int k = 0; k < vector->nnz; ++k)
  {
    SLEQP_CALL(sleqp_sparse_matrix_push(matrix,
                                        vector->indices[k],
                                        column,
                                        vector->data[k] * scale));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE rayleigh_quotient(SleqpSparseVec* direction,
                                       SleqpSparseVec* product,
                                       double eps,
                                       double* rayleigh_factor)
{
  double dir_norm = sleqp_sparse_vector_norm_sq(direction);

  if(sleqp_is_zero(dir_norm, eps))
  {
    return SLEQP_OKAY;
  }

  double dot;

  SLEQP_CALL(sleqp_sparse_vector_dot(direction,
                                     product,
                                     &dot));

  *rayleigh_factor = SLEQP_ABS(dot) / dir_norm;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE trlib_get_status_string(int value,
                                             const char** message)
{
  switch(value)
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

static SLEQP_RETCODE check_optimality(SleqpTRSolver* data,
                                      SleqpAugJacobian* jacobian,
                                      SleqpSparseVec* multipliers,
                                      SleqpSparseVec* gradient,
                                      SleqpSparseVec* newton_step,
                                      bool bdry_sol,
                                      bool* is_optimal)
{
  SleqpProblem* problem = data->problem;
  SleqpFunc* func = problem->func;

  const double eps = sleqp_params_get_eps(data->params);
  const double opt_eps = sleqp_params_get_optimality_tolerance(data->params);
  const double zero_eps = sleqp_params_get_zero_eps(data->params);
  const double one = 1.;

  (*is_optimal) = true;

  SLEQP_CALL(sleqp_func_hess_prod(func,
                                  &one,
                                  newton_step,
                                  multipliers,
                                  data->Hp));

  SLEQP_CALL(sleqp_sparse_vector_add(gradient,
                                     data->Hp,
                                     zero_eps,
                                     data->g));

  SLEQP_CALL(sleqp_aug_jacobian_projection(jacobian,
                                           data->g,
                                           data->v,
                                           NULL));

  const double vnorm = sleqp_sparse_vector_norm(data->v);

  if(bdry_sol)
  {
    if(sleqp_is_zero(vnorm, opt_eps))
    {
      return SLEQP_OKAY;
    }

    double dot;

    SLEQP_CALL(sleqp_sparse_vector_dot(data->v, newton_step, &dot));

    dot *= -1.;

    const double snorm = sleqp_sparse_vector_norm(newton_step);

    assert(snorm > 0.);

    if(!sleqp_is_eq(dot, vnorm*snorm, opt_eps))
    {
      (*is_optimal) = false;
    }

  }
  else
  {
    if(!sleqp_is_zero(vnorm, opt_eps))
    {
      (*is_optimal) = false;
    }

  }

  return SLEQP_OKAY;
}


static SLEQP_RETCODE tr_solve_loop(SleqpTRSolver* data,
                                   SleqpAugJacobian* jacobian,
                                   SleqpSparseVec* multipliers,
                                   SleqpSparseVec* gradient,
                                   double trust_radius,
                                   trlib_int_t* trlib_ret)
{
  SleqpProblem* problem = data->problem;
  SleqpFunc* func = problem->func;

  const double inf = sleqp_infinity();

  const double eps = sleqp_params_get_eps(data->params);
  const double zero_eps = sleqp_params_get_zero_eps(data->params);
  const double rel_tol = sleqp_params_get_newton_relative_tolerance(data->params);

  trlib_int_t equality = 0;
  trlib_int_t maxlanczos = data->trlib_maxiter;
  trlib_int_t ctl_invariant = 0;
  trlib_int_t refine = 1;
  trlib_int_t verbose = (sleqp_log_level() >= SLEQP_LOG_DEBUG);
  trlib_int_t unicode = 0;
  trlib_flt_t tol_rel_i = rel_tol;
  trlib_flt_t tol_abs_i = zero_eps;
  trlib_flt_t tol_rel_b = rel_tol;
  trlib_flt_t tol_abs_b = zero_eps;
  trlib_flt_t obj_lo = -inf;
  trlib_int_t convexify = 1;
  trlib_int_t earlyterm = 1;

  trlib_int_t init = 0;

  trlib_flt_t v_dot_g = 0.0, p_dot_Hp = 0.0, g_dot_g = 0.0, flt1, flt2, flt3;

  trlib_int_t action, ityp;

  //trlib_flt_t zero = TRLIB_EPS*TRLIB_EPS;

  trlib_flt_t zero = zero_eps;

  trlib_int_t iter = 0;

  trlib_int_t *iwork = data->trlib_iwork;
  trlib_flt_t *fwork = data->trlib_fwork;

  trlib_int_t *timing = data->trlib_timinig;

  trlib_int_t maxiter = data->trlib_maxiter;

  SLEQP_CALL(sleqp_timer_start(data->timer));

  {
    init = TRLIB_CLS_INIT;
    trlib_krylov_prepare_memory(maxiter, fwork);
  }

  SLEQP_CALL(sleqp_sparse_matrix_clear(data->Q));

  SLEQP_CALL(sleqp_sparse_vector_clear(data->p));
  SLEQP_CALL(sleqp_sparse_vector_clear(data->Hp));

  char* prefix = "trlib: ";
  FILE* fout = stderr;

  double max_rayleigh = -inf;
  double min_rayleigh = inf;

  while(1)
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

    const double one = 1.;

    init = 0;

    switch(action)
    {
    case TRLIB_CLA_TRIVIAL:
    {
      break;
    }
    case TRLIB_CLA_INIT:
    {
      // reset s to 0
      SLEQP_CALL(sleqp_sparse_vector_clear(data->s));

      SLEQP_CALL(sleqp_sparse_vector_copy(gradient, data->g));

      SLEQP_CALL(sleqp_sparse_vector_clear(data->gm));

      SLEQP_CALL(sleqp_aug_jacobian_projection(jacobian,
                                               data->g,
                                               data->v,
                                               NULL));

      SLEQP_CALL(sleqp_sparse_vector_copy(data->v, data->p));

      {
        double v_norm = sleqp_sparse_vector_norm(data->v);

        // We can stop directly, if the projected gradient is zero
        if(sleqp_is_zero(v_norm, zero_eps))
        {
          return SLEQP_OKAY;
        }

      }

      SLEQP_CALL(sleqp_sparse_vector_scale(data->p, -1.));

      g_dot_g = sleqp_sparse_vector_norm_sq(data->g);

      SLEQP_CALL(sleqp_sparse_vector_dot(data->v,
                                         data->g,
                                         &v_dot_g));

      SLEQP_CALL(sleqp_func_hess_prod(func,
                                      &one,
                                      data->p,
                                      multipliers,
                                      data->Hp));

      if(collect_rayleigh)
      {
        double cur_rayleigh;

        SLEQP_CALL(rayleigh_quotient(data->p, data->Hp, zero_eps, &cur_rayleigh));

        max_rayleigh = SLEQP_MAX(max_rayleigh, cur_rayleigh);
        min_rayleigh = SLEQP_MIN(min_rayleigh, cur_rayleigh);
      }

      SLEQP_CALL(sleqp_sparse_vector_dot(data->p, data->Hp, &p_dot_Hp));

      const int num_rows = sleqp_sparse_matrix_get_num_rows(data->Q);

      SLEQP_CALL(sleqp_sparse_matrix_resize(data->Q, num_rows, 0));
      SLEQP_CALL(sleqp_sparse_matrix_clear(data->Q));

      assert(sleqp_sparse_matrix_get_num_cols(data->Q) == 0);

      //assert(v_dot_g > 0);

      double scale = sqrt(v_dot_g);

      if(!sleqp_is_pos(scale, zero_eps))
      {
        return SLEQP_OKAY;
      }

      scale = 1./ scale;

      SLEQP_CALL(matrix_push_column(data->Q,
                                    data->v,
                                    scale));

      assert(sleqp_sparse_matrix_valid(data->Q));

      break;
    }
    case TRLIB_CLA_RETRANSF:
    {
      SLEQP_CALL(sleqp_sparse_vector_from_raw(data->h,
                                              fwork + data->trlib_h_pointer,
                                              iter + 1,
                                              zero_eps));

      SLEQP_CALL(sleqp_sparse_matrix_vector_product(data->Q, data->h, data->dense_cache));

      SLEQP_CALL(sleqp_sparse_vector_from_raw(data->s,
                                              data->dense_cache,
                                              problem->num_variables,
                                              zero_eps));


      break;
    }
    case TRLIB_CLA_UPDATE_STATIO:
    {
      if(ityp == TRLIB_CLT_CG)
      {
        SLEQP_CALL(sleqp_sparse_vector_add_scaled(data->s,
                                                  data->p,
                                                  1.,
                                                  flt1,
                                                  zero_eps,
                                                  data->sparse_cache));

        SLEQP_CALL(sleqp_sparse_vector_copy(data->sparse_cache, data->s));
      }
      break;
    }
    case TRLIB_CLA_UPDATE_GRAD:
    {
      if(ityp == TRLIB_CLT_CG)
      {
        if(iter + 1 == sleqp_sparse_matrix_get_num_cols(data->Q))
        {
          const int num_rows = sleqp_sparse_matrix_get_num_rows(data->Q);
          SLEQP_CALL(sleqp_sparse_matrix_resize(data->Q, num_rows, iter));
        }

        assert(iter == sleqp_sparse_matrix_get_num_cols(data->Q));

        SLEQP_CALL(matrix_push_column(data->Q,
                                      data->v,
                                      flt2));

        assert(iter + 1 == sleqp_sparse_matrix_get_num_cols(data->Q));
        assert(sleqp_sparse_matrix_valid(data->Q));

        SLEQP_CALL(sleqp_sparse_vector_copy(data->g, data->gm));

        SLEQP_CALL(sleqp_sparse_vector_add_scaled(data->g,
                                                  data->Hp,
                                                  1.,
                                                  flt1,
                                                  zero_eps,
                                                  data->sparse_cache));

        SLEQP_CALL(sleqp_sparse_vector_copy(data->sparse_cache, data->g));

        SLEQP_CALL(sleqp_aug_jacobian_projection(jacobian,
                                                 data->g,
                                                 data->v,
                                                 NULL));

        g_dot_g = sleqp_sparse_vector_norm_sq(data->g);

        SLEQP_CALL(sleqp_sparse_vector_dot(data->v,
                                           data->g,
                                           &v_dot_g));

      }
      if(ityp == TRLIB_CLT_L)
      {
        SLEQP_CALL(sleqp_sparse_vector_add_scaled(data->Hp,
                                                  data->g,
                                                  1.,
                                                  flt1,
                                                  zero_eps,
                                                  data->sparse_cache));

        SLEQP_CALL(sleqp_sparse_vector_add_scaled(data->sparse_cache,
                                                  data->gm,
                                                  1.,
                                                  flt2,
                                                  zero_eps,
                                                  data->s));

        SLEQP_CALL(sleqp_sparse_vector_copy(data->g, data->gm));

        SLEQP_CALL(sleqp_sparse_vector_scale(data->gm, flt3));

        SLEQP_CALL(sleqp_sparse_vector_copy(data->s, data->g));

        SLEQP_CALL(sleqp_aug_jacobian_projection(jacobian,
                                                 data->g,
                                                 data->v,
                                                 NULL));

        SLEQP_CALL(sleqp_sparse_vector_dot(data->v,
                                           data->g,
                                           &v_dot_g));

      }

      break;
    }
    case TRLIB_CLA_CONV_HARD:
    {
      SLEQP_CALL(sleqp_func_hess_prod(func, &one,
                                      data->s,
                                      multipliers,
                                      data->sparse_cache));

      if(collect_rayleigh)
      {
        double cur_rayleigh;

        SLEQP_CALL(rayleigh_quotient(data->s, data->sparse_cache, zero_eps, &cur_rayleigh));

        max_rayleigh = SLEQP_MAX(max_rayleigh, cur_rayleigh);
        min_rayleigh = SLEQP_MIN(min_rayleigh, cur_rayleigh);
      }

      SLEQP_CALL(sleqp_sparse_vector_add(data->sparse_cache, gradient, zero_eps, data->l));

      SLEQP_CALL(sleqp_sparse_vector_add_scaled(data->l,
                                                data->s,
                                                1.,
                                                flt1,
                                                zero_eps,
                                                data->h_lhs));

      SLEQP_CALL(sleqp_aug_jacobian_projection(jacobian,
                                               data->l,
                                               data->sparse_cache,
                                               NULL));

      SLEQP_CALL(sleqp_sparse_vector_add_scaled(data->sparse_cache,
                                                data->s,
                                                1.,
                                                flt1,
                                                zero_eps,
                                                data->h_rhs));

      SLEQP_CALL(sleqp_sparse_vector_dot(data->h_lhs, data->h_rhs, &v_dot_g));

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
      if(ityp == TRLIB_CLT_CG)
      {
        assert(flt1 == -1.);

        SLEQP_CALL(sleqp_sparse_vector_add_scaled(data->v,
                                                  data->p,
                                                  flt1,
                                                  flt2,
                                                  zero_eps,
                                                  data->sparse_cache));

        SLEQP_CALL(sleqp_sparse_vector_copy(data->sparse_cache, data->p));

        SLEQP_CALL(sleqp_func_hess_prod(func, &one,
                                        data->p,
                                        multipliers,
                                        data->Hp));

        if(collect_rayleigh)
        {
          double cur_rayleigh;

          SLEQP_CALL(rayleigh_quotient(data->p, data->Hp, zero_eps, &cur_rayleigh));

          max_rayleigh = SLEQP_MAX(max_rayleigh, cur_rayleigh);
          min_rayleigh = SLEQP_MIN(min_rayleigh, cur_rayleigh);
        }

        SLEQP_CALL(sleqp_sparse_vector_dot(data->p, data->Hp, &p_dot_Hp));

      }
      else if(ityp == TRLIB_CLT_L)
      {
        assert(flt2 == 0.);

        const int num_cols = sleqp_sparse_matrix_get_num_cols(data->Q);

        if(iter + 1 == num_cols)
        {
          const int num_rows = sleqp_sparse_matrix_get_num_rows(data->Q);
          SLEQP_CALL(sleqp_sparse_matrix_resize(data->Q, num_rows, iter));
          assert(sleqp_sparse_matrix_valid(data->Q));
        }

        assert(iter == sleqp_sparse_matrix_get_num_cols(data->Q));

        SLEQP_CALL(sleqp_sparse_vector_add_scaled(data->v,
                                                  data->p,
                                                  flt1,
                                                  flt2,
                                                  zero_eps,
                                                  data->sparse_cache));

        SLEQP_CALL(sleqp_sparse_vector_copy(data->sparse_cache, data->p));

        SLEQP_CALL(sleqp_func_hess_prod(func, &one,
                                        data->p,
                                        multipliers,
                                        data->Hp));

        if(collect_rayleigh)
        {
          double cur_rayleigh;

          SLEQP_CALL(rayleigh_quotient(data->p, data->Hp, zero_eps, &cur_rayleigh));

          max_rayleigh = SLEQP_MAX(max_rayleigh, cur_rayleigh);
          min_rayleigh = SLEQP_MIN(min_rayleigh, cur_rayleigh);
        }

        SLEQP_CALL(sleqp_sparse_vector_dot(data->p, data->Hp, &p_dot_Hp));

        SLEQP_CALL(matrix_push_column(data->Q,
                                      data->p,
                                      1.));

        assert(sleqp_sparse_matrix_valid(data->Q));

      }
      break;
    }
    case TRLIB_CLA_OBJVAL:
    {
      double lin_term, quad_term;

      SLEQP_CALL(sleqp_func_hess_bilinear(func,
                                          &one,
                                          data->s,
                                          multipliers,
                                          &quad_term));

      SLEQP_CALL(sleqp_sparse_vector_dot(data->s,
                                         data->g,
                                         &lin_term));

      g_dot_g = lin_term + .5* quad_term;

      break;
    }
    default:
    {
      sleqp_log_error("Invalid trlib action requested: %d", action);
      return SLEQP_INTERNAL_ERROR;
    }
    }

    if((*trlib_ret) < TRLIB_CLR_CONTINUE)
    {
      break;
    }

    if(data->time_limit != SLEQP_NONE && sleqp_timer_elapsed(data->timer) >= data->time_limit)
    {
      break;
    }
  }

  if(collect_rayleigh && !sleqp_is_zero(min_rayleigh, zero_eps))
  {
    const double cond_bound = SLEQP_ABS(max_rayleigh) / SLEQP_ABS(min_rayleigh);

    sleqp_log_info("Spectrum bound: %f / %f, condition bound: %f",
                   min_rayleigh,
                   max_rayleigh,
                   cond_bound);
  }
}

SLEQP_RETCODE sleqp_tr_solver_solve(SleqpTRSolver* data,
                                    SleqpAugJacobian* jacobian,
                                    SleqpSparseVec* multipliers,
                                    SleqpSparseVec* gradient,
                                    SleqpSparseVec* newton_step,
                                    double trust_radius)
{
  const double eps = sleqp_params_get_eps(data->params);
  const double zero_eps = sleqp_params_get_zero_eps(data->params);

  trlib_int_t ret = 0;

  SLEQP_CALL(sleqp_sparse_vector_clear(newton_step));

  SLEQP_CALL(tr_solve_loop(data,
                           jacobian,
                           multipliers,
                           gradient,
                           trust_radius,
                           &ret));

  // We may loose some orthogonality in the process, reproject to
  // be sure
  SLEQP_CALL(sleqp_aug_jacobian_projection(jacobian,
                                           data->s,
                                           newton_step,
                                           NULL));

  SLEQP_CALL(sleqp_timer_stop(data->timer));

  // check optimality

  const bool solved_optimally = (ret == TRLIB_CLR_CONV_BOUND) ||
    (ret == TRLIB_CLR_CONV_INTERIOR);

  const bool converged_bdry = (ret == TRLIB_CLR_CONV_BOUND);

#if !defined(NDEBUG)
  if(solved_optimally)
  {
    bool tr_subproblem_is_optimal;

    SLEQP_CALL(check_optimality(data,
                                jacobian,
                                multipliers,
                                gradient,
                                newton_step,
                                converged_bdry,
                                &tr_subproblem_is_optimal));

    assert(tr_subproblem_is_optimal);
  }

  // TODO: Choose appropriate tolerance
  // assert(sleqp_is_leq(sleqp_sparse_vector_norm(newton_step), trust_radius, eps));

#endif

  const char* trlib_status_string;

  SLEQP_CALL(trlib_get_status_string(ret, &trlib_status_string));

  if(!solved_optimally)
  {
    sleqp_log_warn("Failed to solve trust region subproblem, reason: %s",
                   trlib_status_string);
  }

  const bool failure_occured =
    (ret == TRLIB_CLR_FAIL_FACTOR) ||
    (ret == TRLIB_CLR_FAIL_LINSOLVE) ||
    (ret == TRLIB_CLR_FAIL_NUMERIC) ||
    (ret == TRLIB_CLR_FAIL_TTR);

  if(failure_occured)
  {
    return SLEQP_INTERNAL_ERROR;
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_tr_solver_capture(SleqpTRSolver* data)
{
  ++data->refcount;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_tr_solver_release(SleqpTRSolver** star)
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
