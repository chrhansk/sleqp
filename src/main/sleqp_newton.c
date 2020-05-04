#include "sleqp_newton.h"

#include <math.h>
#include <trlib.h>

#include "sleqp.h"
#include "sleqp_cmp.h"
#include "sleqp_iterate.h"
#include "sleqp_log.h"
#include "sleqp_mem.h"
#include "sleqp_timer.h"

static const bool collect_rayleigh = false;

struct SleqpNewtonData
{
  SleqpProblem* problem;
  SleqpParams* params;

  // trlib-related data:
  trlib_int_t trlib_maxiter;
  trlib_int_t trlib_h_pointer;
  trlib_int_t *trlib_iwork;
  trlib_flt_t *trlib_fwork;

  SleqpSparseVec* initial_solution;

  SleqpSparseVec* lower_diff;
  SleqpSparseVec* upper_diff;

  SleqpSparseVec* violated_multipliers;
  SleqpSparseVec* multipliers;

  SleqpSparseVec* gradient;

  SleqpSparseVec* initial_rhs;
  SleqpSparseVec* initial_hessian_product;
  SleqpSparseVec* jacobian_product;

  SleqpSparseVec* g;
  SleqpSparseVec* gm;

  SleqpSparseVec* l;
  SleqpSparseVec* h_lhs;
  SleqpSparseVec* h_rhs;

  SleqpSparseVec* h;

  SleqpSparseVec* v;
  SleqpSparseVec* p;
  SleqpSparseVec* Hp;

  double* dense_cache;
  SleqpSparseVec* sparse_cache;

  SleqpSparseMatrix* Q;

  SleqpTimer* timer;
  double time_limit;
};

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

SLEQP_RETCODE sleqp_newton_data_create(SleqpNewtonData** star,
                                       SleqpProblem* problem,
                                       SleqpParams* params)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpNewtonData* data = *star;

  *data = (SleqpNewtonData) {0};

  data->problem = problem;
  data->params = params;

  data->trlib_maxiter = problem->num_variables;

  trlib_int_t iwork_size, fwork_size;


  trlib_krylov_memory_size(data->trlib_maxiter,
                           &iwork_size,
                           &fwork_size,
                           &data->trlib_h_pointer);

  SLEQP_CALL(sleqp_calloc(&data->trlib_iwork, iwork_size));
  SLEQP_CALL(sleqp_calloc(&data->trlib_fwork, fwork_size));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->initial_solution,
                                        problem->num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->lower_diff, problem->num_variables, 0));
  SLEQP_CALL(sleqp_sparse_vector_create(&data->upper_diff, problem->num_variables, 0));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->violated_multipliers, problem->num_constraints, 0));
  SLEQP_CALL(sleqp_sparse_vector_create(&data->multipliers, problem->num_constraints, 0));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->gradient, problem->num_variables, 0));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->initial_rhs, problem->num_constraints, 0));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->initial_hessian_product,
                                        problem->num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->jacobian_product,
                                        problem->num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->g, problem->num_variables, 0));
  SLEQP_CALL(sleqp_sparse_vector_create(&data->gm, problem->num_variables, 0));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->l, problem->num_variables, 0));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->h_lhs, problem->num_variables, 0));
  SLEQP_CALL(sleqp_sparse_vector_create(&data->h_rhs, problem->num_variables, 0));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->h, problem->num_variables, 0));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->v, problem->num_variables, 0));
  SLEQP_CALL(sleqp_sparse_vector_create(&data->p, problem->num_variables, 0));
  SLEQP_CALL(sleqp_sparse_vector_create(&data->Hp, problem->num_variables, 0));

  SLEQP_CALL(sleqp_calloc(&data->dense_cache, SLEQP_MAX(problem->num_variables,
                                                        problem->num_constraints)));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->sparse_cache, problem->num_variables, 0));

  SLEQP_CALL(sleqp_sparse_matrix_create(&data->Q,
                                        problem->num_variables,
                                        data->trlib_maxiter,
                                        0));

  data->time_limit = -1;

  SLEQP_CALL(sleqp_timer_create(&data->timer));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_newton_set_time_limit(SleqpNewtonData* data,
                                          double time_limit)
{
  data->time_limit = time_limit;
  return SLEQP_OKAY;
}

static SLEQP_RETCODE get_initial_rhs(SleqpNewtonData* data,
                                     SleqpIterate* iterate,
                                     SleqpAugJacobian* jacobian)
{
  SleqpProblem* problem = data->problem;

  SleqpSparseVec* initial_rhs = data->initial_rhs;
  SleqpActiveSet* active_set = iterate->active_set;

  const double eps = sleqp_params_get_eps(data->params);

  const int active_set_size = sleqp_active_set_size(iterate->active_set);

  SLEQP_CALL(sleqp_sparse_vector_resize(initial_rhs, active_set_size));

  {
    SLEQP_CALL(sleqp_sparse_vector_clear(initial_rhs));

    SLEQP_CALL(sleqp_sparse_vector_reserve(initial_rhs, active_set_size));

    SLEQP_CALL(sleqp_sparse_vector_resize(initial_rhs, active_set_size));
  }

  // variables
  {
    SleqpSparseVec* values = iterate->primal;
    SleqpSparseVec* var_lb = problem->var_lb;
    SleqpSparseVec* var_ub = problem->var_ub;

    SleqpSparseVec* lower_diff = data->lower_diff;
    SleqpSparseVec* upper_diff = data->upper_diff;

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(values, var_ub, -1., 1., eps, upper_diff));

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(values, var_lb, -1., 1., eps, lower_diff));

    int k_lower = 0, k_upper = 0;

    while(k_lower < lower_diff->nnz || k_upper < upper_diff->nnz)
    {
      bool valid_lower = (k_lower < lower_diff->nnz);
      bool valid_upper = (k_upper < upper_diff->nnz);

      int i_lower = valid_lower ? lower_diff->indices[k_lower] : lower_diff->dim + 1;
      int i_upper = valid_upper ? upper_diff->indices[k_upper] : upper_diff->dim + 1;

      int i_combined = SLEQP_MIN(i_lower, i_upper);

      valid_lower = valid_lower && (i_lower == i_combined);
      valid_upper = valid_upper && (i_upper == i_combined);

      double lower_value = valid_lower ? lower_diff->data[k_lower] : 0.;
      double upper_value = valid_upper ? upper_diff->data[k_upper] : 0.;

      int i_set = sleqp_active_set_get_variable_index(active_set, i_combined);

      SLEQP_ACTIVE_STATE var_state = sleqp_active_set_get_variable_state(active_set,
                                                                         i_combined);

      if(var_state == SLEQP_ACTIVE_UPPER)
      {
        SLEQP_CALL(sleqp_sparse_vector_push(initial_rhs,
                                            i_set,
                                            upper_value));
      }
      else if(var_state == SLEQP_ACTIVE_LOWER)
      {
        SLEQP_CALL(sleqp_sparse_vector_push(initial_rhs,
                                            i_set,
                                            lower_value));
      }
      else if(var_state == SLEQP_ACTIVE_BOTH)
      {
        assert(sleqp_eq(lower_value, upper_value, eps));

        SLEQP_CALL(sleqp_sparse_vector_push(initial_rhs,
                                            i_set,
                                            lower_value));
      }

      if(i_lower == i_combined)
      {
        ++k_lower;
      }

      if(i_upper == i_combined)
      {
        ++k_upper;
      }

    }
  }

  // constraints
  {
    SleqpSparseVec* values = iterate->cons_val;
    SleqpSparseVec* cons_lb = problem->cons_lb;
    SleqpSparseVec* cons_ub = problem->cons_ub;

    SleqpSparseVec* lower_diff = data->lower_diff;
    SleqpSparseVec* upper_diff = data->upper_diff;

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(values, cons_ub, -1., 1., eps, upper_diff));

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(values, cons_lb, -1., 1., eps, lower_diff));

    int k_lower = 0, k_upper = 0;

    while(k_lower < lower_diff->nnz || k_upper < upper_diff->nnz)
    {
      bool valid_lower = (k_lower < lower_diff->nnz);
      bool valid_upper = (k_upper < upper_diff->nnz);

      int i_lower = valid_lower ? lower_diff->indices[k_lower] : lower_diff->dim + 1;
      int i_upper = valid_upper ? upper_diff->indices[k_upper] : upper_diff->dim + 1;

      int i_combined = SLEQP_MIN(i_lower, i_upper);

      valid_lower = valid_lower && (i_lower == i_combined);
      valid_upper = valid_upper && (i_upper == i_combined);

      double lower_value = valid_lower ? lower_diff->data[k_lower] : 0.;
      double upper_value = valid_upper ? upper_diff->data[k_upper] : 0.;

      int i_set = sleqp_active_set_get_constraint_index(active_set, i_combined);

      SLEQP_ACTIVE_STATE cons_state = sleqp_active_set_get_constraint_state(active_set,
                                                                            i_combined);

      if(cons_state == SLEQP_ACTIVE_UPPER)
      {
        SLEQP_CALL(sleqp_sparse_vector_push(initial_rhs,
                                            i_set,
                                            upper_value));
      }
      else if(cons_state == SLEQP_ACTIVE_LOWER)
      {
        SLEQP_CALL(sleqp_sparse_vector_push(initial_rhs,
                                            i_set,
                                            lower_value));
      }
      else if(cons_state == SLEQP_ACTIVE_BOTH)
      {
        assert(sleqp_eq(lower_value, upper_value, eps));

        SLEQP_CALL(sleqp_sparse_vector_push(initial_rhs,
                                            i_set,
                                            lower_value));
      }

      if(valid_lower)
      {
        ++k_lower;
      }

      if(valid_upper)
      {
        ++k_upper;
      }
    }
  }

  return SLEQP_OKAY;
}


static SLEQP_RETCODE matrix_push_column(SleqpSparseMatrix* matrix,
                                        SleqpSparseVec* vector,
                                        double scale)
{
  SLEQP_CALL(sleqp_sparse_matrix_reserve(matrix, matrix->nnz + vector->nnz));

  const int column = matrix->num_cols;

  ++matrix->num_cols;

  SLEQP_CALL(sleqp_sparse_matrix_push_column(matrix, column));

  for(int k = 0; k < vector->nnz; ++k)
  {
    SLEQP_CALL(sleqp_sparse_matrix_push(matrix,
                                        vector->indices[k],
                                        column,
                                        vector->data[k] * scale));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE matrix_pop_column(SleqpSparseMatrix* matrix)
{
  assert(matrix->num_cols > 0);

  int col_nnz = matrix->cols[matrix->num_cols] - matrix->cols[matrix->num_cols - 1];

  matrix->nnz -= col_nnz;
  --matrix->num_cols;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE rayleigh_quotient(SleqpSparseVec* direction,
                                       SleqpSparseVec* product,
                                       double eps,
                                       double* rayleigh_factor)
{
  double dir_norm = sleqp_sparse_vector_normsq(direction);

  if(sleqp_zero(dir_norm, eps))
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

static SLEQP_RETCODE solve_trust_region_subproblem(SleqpNewtonData* data,
                                                   SleqpIterate* iterate,
                                                   SleqpAugJacobian* jacobian,
                                                   SleqpSparseVec* newton_step,
                                                   double trust_radius)
{
  SleqpProblem* problem = data->problem;
  SleqpFunc* func = problem->func;

  trlib_int_t equality = 0;
  trlib_int_t maxlanczos = 100;
  trlib_int_t ctl_invariant = 0;
  trlib_int_t refine = 1;
  trlib_int_t verbose = (sleqp_log_level() >= SLEQP_LOG_DEBUG);
  trlib_int_t unicode = 0;
  trlib_flt_t tol_rel_i = -2.0;
  trlib_flt_t tol_abs_i = 0.0;
  trlib_flt_t tol_rel_b = -3.0;
  trlib_flt_t tol_abs_b = 0.0;
  trlib_flt_t obj_lo = -1e20;
  trlib_int_t convexify = 1;
  trlib_int_t earlyterm = 1;

  trlib_int_t init = 0;
  trlib_int_t ret = 0;

  trlib_flt_t v_dot_g = 0.0, p_dot_Hp = 0.0, g_dot_g = 0.0, flt1, flt2, flt3;

  trlib_int_t action, ityp;

  //trlib_flt_t zero = TRLIB_EPS*TRLIB_EPS;

  const double zero_eps = sleqp_params_get_zero_eps(data->params);

  trlib_flt_t zero = zero_eps;

  trlib_int_t iter = 0;

  trlib_int_t *iwork = data->trlib_iwork;
  trlib_flt_t *fwork = data->trlib_fwork;

  trlib_int_t maxiter = data->trlib_maxiter;
  //trlib_int_t timing = 0;

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

  double max_rayleigh = -sleqp_infinity();
  double min_rayleigh = sleqp_infinity();

  while(1)
  {
    ret = trlib_krylov_min(init,          // trlib_int_t init
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
                           NULL,          // trlib_int_t *timing
                           //&timing,       // trlib_int_t *timing
                           &action,       // trlib_int_t *action
                           &iter,         // trlib_int_t *iter
                           &ityp,         // trlib_int_t *ityp
                           &flt1,         // trlib_flt_t *flt1
                           &flt2,         // trlib_flt_t *flt2
                           &flt3);        // trlib_flt_t *flt3

    double one = 1.;

    init = 0;

    switch(action)
    {
    case TRLIB_CLS_INIT:
    {
      // reset s to 0
      SLEQP_CALL(sleqp_sparse_vector_clear(newton_step));

      SLEQP_CALL(sleqp_sparse_vector_copy(data->gradient, data->g));

      SLEQP_CALL(sleqp_sparse_vector_clear(data->gm));

      SLEQP_CALL(sleqp_aug_jacobian_projection(jacobian,
                                               data->g,
                                               data->v,
                                               NULL));

      SLEQP_CALL(sleqp_sparse_vector_copy(data->v, data->p));

      {
        double v_norm = sleqp_sparse_vector_norm(data->v);

        // We can stop directly, if the projected gradient is zero
        if(sleqp_zero(v_norm, zero_eps))
        {
          return SLEQP_OKAY;
        }

      }

      SLEQP_CALL(sleqp_sparse_vector_scale(data->p, -1.));

      g_dot_g = sleqp_sparse_vector_normsq(data->g);

      SLEQP_CALL(sleqp_sparse_vector_dot(data->v,
                                         data->g,
                                         &v_dot_g));

      SLEQP_CALL(sleqp_func_hess_prod(func,
                                      &one,
                                      data->p,
                                      data->multipliers,
                                      data->Hp));

      if(collect_rayleigh)
      {
        double cur_rayleigh;

        SLEQP_CALL(rayleigh_quotient(data->p, data->Hp, zero_eps, &cur_rayleigh));

        max_rayleigh = SLEQP_MAX(max_rayleigh, cur_rayleigh);
        min_rayleigh = SLEQP_MIN(min_rayleigh, cur_rayleigh);
      }

      SLEQP_CALL(sleqp_sparse_vector_dot(data->p, data->Hp, &p_dot_Hp));

      SLEQP_CALL(sleqp_sparse_matrix_clear(data->Q));

      data->Q->num_cols = 0;

      //assert(v_dot_g > 0);

      double scale = sqrt(v_dot_g);

      if(!sleqp_pos(scale, zero_eps))
      {
        return SLEQP_OKAY;
      }

      scale = 1./ scale;

      SLEQP_CALL(matrix_push_column(data->Q,
                                    data->v,
                                    scale));

      break;
    }
    case TRLIB_CLA_RETRANSF:
    {
      SLEQP_CALL(sleqp_sparse_vector_from_raw(data->h,
                                              fwork + data->trlib_h_pointer,
                                              iter + 1,
                                              zero_eps));

      SLEQP_CALL(sleqp_sparse_matrix_vector_product(data->Q, data->h, data->dense_cache));

      SLEQP_CALL(sleqp_sparse_vector_from_raw(newton_step,
                                              data->dense_cache,
                                              problem->num_variables,
                                              zero_eps));


      break;
    }
    case TRLIB_CLA_UPDATE_STATIO:
    {
      if(ityp == TRLIB_CLT_CG)
      {
        SLEQP_CALL(sleqp_sparse_vector_add_scaled(newton_step,
                                                  data->p,
                                                  1.,
                                                  flt1,
                                                  zero_eps,
                                                  data->sparse_cache));

        SLEQP_CALL(sleqp_sparse_vector_copy(data->sparse_cache, newton_step));
      }
      break;
    }
    case TRLIB_CLA_UPDATE_GRAD:
    {
      if(ityp == TRLIB_CLT_CG)
      {
        if(iter + 1 == data->Q->num_cols)
        {
          matrix_pop_column(data->Q);
        }

        assert(iter == data->Q->num_cols);

        SLEQP_CALL(matrix_push_column(data->Q,
                                      data->v,
                                      flt2));

        assert(iter + 1 == data->Q->num_cols);

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

        g_dot_g = sleqp_sparse_vector_normsq(data->g);

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
                                                  newton_step));

        SLEQP_CALL(sleqp_sparse_vector_copy(data->g, data->gm));

        SLEQP_CALL(sleqp_sparse_vector_scale(data->gm, flt3));

        SLEQP_CALL(sleqp_sparse_vector_copy(newton_step, data->g));

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
                                      newton_step,
                                      data->multipliers,
                                      data->sparse_cache));

      if(collect_rayleigh)
      {
        double cur_rayleigh;

        SLEQP_CALL(rayleigh_quotient(newton_step, data->sparse_cache, zero_eps, &cur_rayleigh));

        max_rayleigh = SLEQP_MAX(max_rayleigh, cur_rayleigh);
        min_rayleigh = SLEQP_MIN(min_rayleigh, cur_rayleigh);
      }

      SLEQP_CALL(sleqp_sparse_vector_add(data->sparse_cache, data->gradient, zero_eps, data->l));

      SLEQP_CALL(sleqp_sparse_vector_add_scaled(data->l,
                                                newton_step,
                                                1.,
                                                flt1,
                                                zero_eps,
                                                data->h_lhs));

      SLEQP_CALL(sleqp_aug_jacobian_projection(jacobian,
                                               data->l,
                                               data->sparse_cache,
                                               NULL));

      SLEQP_CALL(sleqp_sparse_vector_add_scaled(data->sparse_cache,
                                                newton_step,
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
                                        data->multipliers,
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
      if(ityp == TRLIB_CLT_L)
      {
        assert(flt2 == 0.);

        if(iter + 1 == data->Q->num_cols)
        {
          matrix_pop_column(data->Q);
        }

        assert(iter == data->Q->num_cols);

        SLEQP_CALL(sleqp_sparse_vector_add_scaled(data->v,
                                                  data->p,
                                                  flt1,
                                                  flt2,
                                                  zero_eps,
                                                  data->sparse_cache));

        SLEQP_CALL(sleqp_sparse_vector_copy(data->sparse_cache, data->p));

        SLEQP_CALL(sleqp_func_hess_prod(func, &one,
                                        data->p,
                                        data->multipliers,
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

      }
      break;
    }
    case TRLIB_CLA_OBJVAL:
    {
      double lin_term, quad_term;

      SLEQP_CALL(sleqp_func_hess_bilinear(func,
                                          &one,
                                          newton_step,
                                          data->multipliers,
                                          &quad_term));

      SLEQP_CALL(sleqp_sparse_vector_dot(newton_step,
                                         data->g,
                                         &lin_term));

      g_dot_g = lin_term + .5* quad_term;

      break;
    }
    }

    if(ret < TRLIB_CLR_CONTINUE)
    {
      break;
    }

    if(data->time_limit != -1 && sleqp_timer_elapsed(data->timer) >= data->time_limit)
    {
      break;
    }
  }

  SLEQP_CALL(sleqp_timer_stop(data->timer));

  assert(iter + 1 == data->Q->num_cols);

  const char* trlib_status_string;

  SLEQP_CALL(trlib_get_status_string(ret, &trlib_status_string));

  if(collect_rayleigh && !sleqp_zero(min_rayleigh, zero_eps))
  {
    const double cond_bound = SLEQP_ABS(max_rayleigh) / SLEQP_ABS(min_rayleigh);

    sleqp_log_info("Spectrum bound: %f / %f, condition bound: %f",
                   min_rayleigh,
                   max_rayleigh,
                   cond_bound);
  }

  if(ret < 0)
  {
    if(ret == TRLIB_CLR_PCINDEF)
    {
      sleqp_log_warn("trlib reported indefinite preconditioner");
    }
    else
    {
      sleqp_log_warn("Caught trlib error <%d> (%s)",
                     ret,
                     trlib_status_string);

      //return SLEQP_INTERNAL_ERROR;
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_newton_compute_step(SleqpNewtonData* data,
                                        SleqpIterate* iterate,
                                        SleqpAugJacobian* jacobian,
                                        double trust_radius,
                                        double penalty_parameter,
                                        SleqpSparseVec* newton_step)
{
  SleqpProblem* problem = data->problem;

  const double eps = sleqp_params_get_eps(data->params);
  const double zero_eps = sleqp_params_get_zero_eps(data->params);

  SLEQP_CALL(get_initial_rhs(data, iterate, jacobian));

  SLEQP_CALL(sleqp_aug_jacobian_min_norm_solution(jacobian,
                                                  data->initial_rhs,
                                                  data->initial_solution));

  const double norm_ratio = .8;

  // rescale min norm solution if required.
  {
    double unscaled_norm = sleqp_sparse_vector_norm(data->initial_solution);

    double initial_norm = unscaled_norm;

    double alpha = 1.;

    if(initial_norm != 0.)
    {
      assert(initial_norm > 0.);

      alpha = (norm_ratio * trust_radius) / (initial_norm);

      alpha = SLEQP_MIN(alpha, 1.);


      if(sleqp_eq(alpha, 1., eps))
      {
        // no scaling required...

        const double initial_norm_sq = initial_norm * initial_norm;

        const double trust_radius_sq = trust_radius * trust_radius;

        //assert(sleqp_lt(initial_norm_sq, trust_radius_sq, eps));

        trust_radius = sqrt(trust_radius_sq - initial_norm_sq);
      }
      else
      {

        SLEQP_CALL(sleqp_sparse_vector_scale(data->initial_solution, alpha));

        // we know that the scaled initial solution
        // has norm equal to norm_ratio * trust_radius
        trust_radius *= sqrt(1. - norm_ratio * norm_ratio);
      }
    }
  }

  int active_set_size = sleqp_active_set_size(iterate->active_set);

  // in either case the only feasible solution is the zero vector
  if(sleqp_zero(trust_radius, zero_eps) || problem->num_variables == active_set_size)
  {
    SLEQP_CALL(sleqp_sparse_vector_copy(data->initial_solution, newton_step));

    return SLEQP_OKAY;
  }

  // compute the EQP multipliers. Given as the sum of the
  // estimated multipliers and penalties on violated
  // constraints
  {
    SLEQP_CALL(sleqp_get_violated_multipliers(problem,
                                              iterate->primal,
                                              iterate->cons_val,
                                              penalty_parameter,
                                              data->violated_multipliers,
                                              iterate->active_set,
                                              zero_eps));

    SLEQP_CALL(sleqp_sparse_vector_add(data->violated_multipliers,
                                       iterate->cons_dual,
                                       zero_eps,
                                       data->multipliers));
  }

  {
    double one = 1.;

    SLEQP_CALL(sleqp_func_hess_prod(problem->func,
                                    &one,
                                    data->initial_solution,
                                    data->multipliers,
                                    data->initial_hessian_product));

    SLEQP_CALL(sleqp_sparse_vector_add(iterate->func_grad,
                                       iterate->vars_dual,
                                       zero_eps,
                                       data->gradient));

    SLEQP_CALL(sleqp_sparse_vector_add(data->gradient,
                                       data->initial_hessian_product,
                                       zero_eps,
                                       data->sparse_cache));

    SLEQP_CALL(sleqp_sparse_matrix_trans_vector_product(iterate->cons_jac,
                                                        data->multipliers,
                                                        zero_eps,
                                                        data->jacobian_product));

    SLEQP_CALL(sleqp_sparse_vector_add(data->sparse_cache,
                                       data->jacobian_product,
                                       zero_eps,
                                       data->gradient));
  }

  SLEQP_CALL(solve_trust_region_subproblem(data,
                                           iterate,
                                           jacobian,
                                           newton_step,
                                           trust_radius));

  SLEQP_CALL(sleqp_sparse_vector_copy(newton_step, data->sparse_cache));

  SLEQP_CALL(sleqp_sparse_vector_add(data->sparse_cache,
                                     data->initial_solution,
                                     zero_eps,
                                     newton_step));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_newton_data_free(SleqpNewtonData** star)
{
  SleqpNewtonData* data = *star;

  if(!data)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_timer_free(&data->timer));

  SLEQP_CALL(sleqp_sparse_matrix_free(&data->Q));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->sparse_cache));
  sleqp_free(&data->dense_cache);

  SLEQP_CALL(sleqp_sparse_vector_free(&data->Hp));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->p));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->v));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->h));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->h_rhs));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->h_lhs));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->l));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->gm));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->g));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->jacobian_product));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->initial_hessian_product));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->initial_rhs));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->gradient));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->multipliers));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->violated_multipliers));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->upper_diff));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->lower_diff));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->initial_solution));

  sleqp_free(&data->trlib_fwork);
  sleqp_free(&data->trlib_iwork);

  sleqp_free(star);

  return SLEQP_OKAY;
}
