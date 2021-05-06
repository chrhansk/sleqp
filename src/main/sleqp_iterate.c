#include "sleqp_iterate.h"

#include <math.h>

#include "sleqp_assert.h"
#include "sleqp_cmp.h"
#include "sleqp_feas.h"
#include "sleqp_mem.h"

struct SleqpIterate
{
  int refcount;
  /**
   * The current point. Has dimension = num_variables.
   **/
  SleqpSparseVec* primal;

  /**
   * The current function value
   **/
  double func_val;

  /**
   * The current function gradient. Has dimension = num_variables.
   **/
  SleqpSparseVec* func_grad;

  /**
   * The current constraint values. Has dimension = num_constraints.
   **/
  SleqpSparseVec* cons_val;

  /**
   * The Jacobian of the constraitns at the current iterate.
   * Has num_constraints many rows, num_variables many columns.
   */
  SleqpSparseMatrix* cons_jac;

  /**
   * The current working set.
   **/
  SleqpWorkingSet* working_set;

  /**
   * The dual values of the constraints. Has dimension = num_constraints.
   */
  SleqpSparseVec* cons_dual;

  /**
   * The dual values of the variable bounds. Has dimension = num_variables.
   */
  SleqpSparseVec* vars_dual;

};

SLEQP_RETCODE sleqp_iterate_create(SleqpIterate** star,
                                   SleqpProblem* problem,
                                   SleqpSparseVec* x)
{
  SLEQP_CALL(sleqp_malloc(star));

  sleqp_assert_msg(sleqp_sparse_vector_is_valid(x), "Invalid primal values");

  sleqp_assert_msg(sleqp_sparse_vector_is_finite(x), "Infinite primal values");

  SleqpIterate* iterate = *star;

  *iterate = (SleqpIterate){0};

  iterate->refcount = 1;

  int num_variables = problem->num_variables;
  int num_constraints = problem->num_constraints;

  SLEQP_CALL(sleqp_sparse_vector_create(&iterate->primal,
                                        num_variables,
                                        x->nnz));

  SLEQP_CALL(sleqp_sparse_vector_copy(x, iterate->primal));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&iterate->func_grad,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&iterate->cons_val,
                                              num_constraints));

  SLEQP_CALL(sleqp_sparse_matrix_create(&iterate->cons_jac,
                                        num_constraints,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_working_set_create(&iterate->working_set,
                                      problem));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&iterate->cons_dual,
                                              num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&iterate->vars_dual,
                                              num_variables));

  return SLEQP_OKAY;
}


SleqpSparseVec* sleqp_iterate_get_primal(SleqpIterate* iterate)
{
  return iterate->primal;
}


double sleqp_iterate_get_func_val(SleqpIterate* iterate)
{
  return iterate->func_val;
}

SLEQP_RETCODE sleqp_iterate_set_func_val(SleqpIterate* iterate,
                                         double value)
{
  iterate->func_val = value;
  return SLEQP_OKAY;
}


SleqpSparseVec* sleqp_iterate_get_func_grad(SleqpIterate* iterate)
{
  return iterate->func_grad;
}


SleqpSparseVec* sleqp_iterate_get_cons_val(SleqpIterate* iterate)
{
  return iterate->cons_val;
}

SleqpSparseMatrix* sleqp_iterate_get_cons_jac(SleqpIterate* iterate)
{
  return iterate->cons_jac;
}

SleqpWorkingSet* sleqp_iterate_get_working_set(SleqpIterate* iterate)
{
  return iterate->working_set;
}

SleqpSparseVec* sleqp_iterate_get_cons_dual(SleqpIterate* iterate)
{
  return iterate->cons_dual;
}

SleqpSparseVec* sleqp_iterate_get_vars_dual(SleqpIterate* iterate)
{
  return iterate->vars_dual;
}

static
SLEQP_RETCODE slack_residuum(const SleqpSparseVec* v,
                             const SleqpSparseVec* lb,
                             const SleqpSparseVec* ub,
                             const SleqpSparseVec* d,
                             double* residuum)
{
  const int dim = v->dim;

  *residuum = 0.;

  {
    {
      int k_v = 0, k_ub = 0, k_lb = 0, k_d = 0;

      int dim = v->dim;

      while((k_v < v->nnz || k_lb < lb->nnz || k_ub < ub->nnz) && k_d < d->nnz)
      {
        bool valid_v = k_v < v->nnz;
        bool valid_ub = k_ub < ub->nnz;
        bool valid_lb = k_lb < ub->nnz;
        bool valid_d = k_d < d->nnz;

        const int i_v = valid_v ? v->indices[k_v] : dim + 1;
        const int i_lb = valid_lb ? lb->indices[k_lb] : dim + 1;
        const int i_ub = valid_ub ? ub->indices[k_ub] : dim + 1;
        const int i_d = valid_d ? d->indices[k_d] : dim + 1;

        int i_combined = SLEQP_MIN(i_lb, i_ub);
        i_combined = SLEQP_MIN(i_combined, i_d);
        i_combined = SLEQP_MIN(i_combined, i_v);

        valid_v = valid_v && (i_v == i_combined);
        valid_lb = valid_lb && (i_lb == i_combined);
        valid_ub = valid_ub && (i_ub == i_combined);
        valid_d = valid_d && (i_d == i_combined);

        const double v_val = valid_v ? v->data[k_v] : 0.;
        const double lb_val = valid_lb ? lb->data[k_lb] : 0.;
        const double ub_val = valid_ub ? ub->data[k_ub] : 0.;
        const double d_val = valid_d ? d->data[k_d] : 0.;

        double current_residuum = 0.;

        // use signs of dual variables with respect to working set
        if(d_val >= 0.)
        {
          current_residuum = SLEQP_MAX(ub_val - v_val, 0.) * d_val;
        }
        else
        {
          current_residuum = SLEQP_MAX(v_val - lb_val, 0.) * d_val;
        }

        *residuum = SLEQP_MAX(*residuum, SLEQP_ABS(current_residuum));

        if(valid_v)
        {
          ++k_v;
        }

        if(valid_lb)
        {
          ++k_lb;
        }

        if(valid_ub)
        {
          ++k_ub;
        }

        if(valid_d)
        {
          ++k_d;
        }
      }
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_iterate_slackness_residuum(SleqpProblem* problem,
                                               SleqpIterate* iterate,
                                               double* slackness_residuum)
{
  (*slackness_residuum) = 0.;

  double current_residuum;

  SLEQP_CALL(slack_residuum(sleqp_iterate_get_cons_val(iterate),
                            problem->var_lb,
                            problem->var_ub,
                            sleqp_iterate_get_vars_dual(iterate),
                            &current_residuum));

  (*slackness_residuum) = SLEQP_MAX((*slackness_residuum), current_residuum);

  SLEQP_CALL(slack_residuum(sleqp_iterate_get_primal(iterate),
                            problem->cons_lb,
                            problem->cons_ub,
                            sleqp_iterate_get_cons_dual(iterate),
                            &current_residuum));

  (*slackness_residuum) = SLEQP_MAX((*slackness_residuum), current_residuum);


  return SLEQP_OKAY;
}

static
SLEQP_RETCODE slack_residuals(const SleqpSparseVec* v,
                              const SleqpSparseVec* lb,
                              const SleqpSparseVec* ub,
                              const SleqpSparseVec* d,
                              SleqpSparseVec* r,
                              double zero_eps)
{
  const int dim = v->dim;

  SLEQP_CALL(sleqp_sparse_vector_clear(r));

  SLEQP_CALL(sleqp_sparse_vector_resize(r, dim));
  SLEQP_CALL(sleqp_sparse_vector_reserve(r, dim));

  {
    {
      int k_v = 0, k_ub = 0, k_lb = 0, k_d = 0;

      int dim = v->dim;

      while((k_v < v->nnz || k_lb < lb->nnz || k_ub < ub->nnz) && k_d < d->nnz)
      {
        bool valid_v = k_v < v->nnz;
        bool valid_ub = k_ub < ub->nnz;
        bool valid_lb = k_lb < ub->nnz;
        bool valid_d = k_d < d->nnz;

        const int i_v = valid_v ? v->indices[k_v] : dim + 1;
        const int i_lb = valid_lb ? lb->indices[k_lb] : dim + 1;
        const int i_ub = valid_ub ? ub->indices[k_ub] : dim + 1;
        const int i_d = valid_d ? d->indices[k_d] : dim + 1;

        int i_combined = SLEQP_MIN(i_lb, i_ub);
        i_combined = SLEQP_MIN(i_combined, i_d);
        i_combined = SLEQP_MIN(i_combined, i_v);

        valid_v = valid_v && (i_v == i_combined);
        valid_lb = valid_lb && (i_lb == i_combined);
        valid_ub = valid_ub && (i_ub == i_combined);
        valid_d = valid_d && (i_d == i_combined);

        const double v_val = valid_v ? v->data[k_v] : 0.;
        const double lb_val = valid_lb ? lb->data[k_lb] : 0.;
        const double ub_val = valid_ub ? ub->data[k_ub] : 0.;
        const double d_val = valid_d ? d->data[k_d] : 0.;

        double current_residuum = 0.;

        // use signs of dual variables with respect to working set
        if(d_val >= 0.)
        {
          current_residuum = SLEQP_MAX(ub_val - v_val, 0.) * d_val;
        }
        else
        {
          current_residuum = SLEQP_MAX(v_val - lb_val, 0.) * d_val;
        }

        if(!sleqp_is_zero(current_residuum, zero_eps))
        {
          SLEQP_CALL(sleqp_sparse_vector_push(r,
                                              i_combined,
                                              current_residuum));
        }

        if(valid_v)
        {
          ++k_v;
        }

        if(valid_lb)
        {
          ++k_lb;
        }

        if(valid_ub)
        {
          ++k_ub;
        }

        if(valid_d)
        {
          ++k_d;
        }
      }
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_iterate_vars_slackness_residuals(SleqpProblem* problem,
                                                     SleqpIterate* iterate,
                                                     SleqpSparseVec* residuals,
                                                     double zero_eps)
{
  const int num_variables = problem->num_variables;

  SLEQP_CALL(slack_residuals(sleqp_iterate_get_cons_val(iterate),
                             problem->var_lb,
                             problem->var_ub,
                             sleqp_iterate_get_vars_dual(iterate),
                             residuals,
                             zero_eps));

  return SLEQP_OKAY;
}


SLEQP_RETCODE sleqp_iterate_cons_slackness_residuals(SleqpProblem* problem,
                                                     SleqpIterate* iterate,
                                                     SleqpSparseVec* residuals,
                                                     double zero_eps)
{
  const int num_constraints = problem->num_constraints;

  SLEQP_CALL(slack_residuals(sleqp_iterate_get_primal(iterate),
                             problem->cons_lb,
                             problem->cons_ub,
                             sleqp_iterate_get_cons_dual(iterate),
                             residuals,
                             zero_eps));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_iterate_feasibility_residuum(SleqpProblem* problem,
                                                 SleqpIterate* iterate,
                                                 double feas_eps,
                                                 double* feasibility_residuum)
{
  SLEQP_CALL(sleqp_violation_inf_norm(problem,
                                      iterate->cons_val,
                                      feas_eps,
                                      feasibility_residuum));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_iterate_get_violated_constraints(SleqpProblem* problem,
                                                     SleqpIterate* iterate,
                                                     int* violated_constraints,
                                                     int* num_violated_constraints,
                                                     double feas_eps)
{
  SLEQP_CALL(sleqp_violated_constraints(problem,
                                        iterate->cons_val,
                                        violated_constraints,
                                        num_violated_constraints,
                                        feas_eps));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
write_stationarity_resiudals_to_cache(SleqpIterate* iterate,
                                      SleqpProblem* problem,
                                      double* cache)
{
  const int num_variables = problem->num_variables;
  const int num_constraints = problem->num_constraints;

  SleqpSparseMatrix* cons_jac = sleqp_iterate_get_cons_jac(iterate);
  SleqpSparseVec* cons_dual = sleqp_iterate_get_cons_dual(iterate);

  SleqpSparseVec* vars_dual = sleqp_iterate_get_vars_dual(iterate);
  SleqpSparseVec* func_grad = sleqp_iterate_get_func_grad(iterate);

  const int num_rows = sleqp_sparse_matrix_get_num_rows(cons_jac);
  const int num_cols = sleqp_sparse_matrix_get_num_cols(cons_jac);

  int* cons_jac_cols = sleqp_sparse_matrix_get_cols(cons_jac);
  int* cons_jac_rows = sleqp_sparse_matrix_get_rows(cons_jac);
  double* cons_jac_data = sleqp_sparse_matrix_get_data(cons_jac);

  assert(num_variables == num_cols);
  assert(num_constraints == num_rows);

  for(int j = 0; j < num_variables; ++j)
  {
    int k_d = 0, k_j = cons_jac_cols[j];

    double sum = 0.;

    while(k_d < cons_dual->nnz && k_j < cons_jac_cols[j + 1])
    {
      int d_idx = cons_dual->indices[k_d];
      int j_idx = cons_jac_rows[k_j];

      if(d_idx < j_idx)
      {
        ++k_d;
      }
      else if(d_idx > j_idx)
      {
        ++k_j;
      }
      else
      {
        sum += cons_dual->data[k_d++] * cons_jac_data[k_j++];
      }
    }

    cache[j] = sum;
  }

  for(int k = 0; k < vars_dual->nnz; ++k)
  {
    cache[vars_dual->indices[k]] += vars_dual->data[k];
  }

  for(int k = 0; k < func_grad->nnz; ++k)
  {
    cache[func_grad->indices[k]] += func_grad->data[k];
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_iterate_stationarity_residuals(SleqpProblem* problem,
                                                   SleqpIterate* iterate,
                                                   double* cache,
                                                   SleqpSparseVec* residuals,
                                                   double zero_eps)
{
  const int num_variables = problem->num_variables;

  SLEQP_CALL(write_stationarity_resiudals_to_cache(iterate, problem, cache));

  SLEQP_CALL(sleqp_sparse_vector_from_raw(residuals,
                                          cache,
                                          num_variables,
                                          zero_eps));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_iterate_stationarity_residuum(SleqpProblem* problem,
                                                  SleqpIterate* iterate,
                                                  double* cache,
                                                  double* stationarity_residuum)
{
  const int num_variables = problem->num_variables;

  SLEQP_CALL(write_stationarity_resiudals_to_cache(iterate, problem, cache));

  (*stationarity_residuum) = 0.;

  for(int j = 0; j < num_variables; ++j)
  {
    (*stationarity_residuum) = SLEQP_MAX((*stationarity_residuum), SLEQP_ABS(cache[j]));
  }

  return SLEQP_OKAY;
}

bool sleqp_iterate_is_feasible(SleqpIterate* iterate,
                               double feasibility_residuum,
                               double feasibility_tolerance)
{
  return feasibility_residuum <= feasibility_tolerance;
}

bool sleqp_iterate_is_optimal(SleqpIterate* iterate,
                              SleqpParams* params,
                              double feasibility_residuum,
                              double slackness_residuum,
                              double stationarity_residuum)
{
  const double feas_eps = sleqp_params_get(params,
                                           SLEQP_PARAM_FEASIBILITY_TOL);

  const double slack_eps = sleqp_params_get(params,
                                            SLEQP_PARAM_SLACKNESS_TOL);

  const double stat_eps = sleqp_params_get(params,
                                           SLEQP_PARAM_STATIONARITY_TOL);

  if(!sleqp_iterate_is_feasible(iterate, feasibility_residuum, feas_eps))
  {
    sleqp_log_debug("Iterate is not feasible, residuum: %e",
                    stationarity_residuum);

    return false;
  }

  double multiplier_norm = 0.;

  {
    multiplier_norm += sleqp_sparse_vector_norm_sq(iterate->cons_dual);
    multiplier_norm += sleqp_sparse_vector_norm_sq(iterate->vars_dual);

    multiplier_norm = sqrt(multiplier_norm);
  }

  if(stationarity_residuum >= stat_eps * (1. + multiplier_norm))
  {
    sleqp_log_debug("Iterate is not optimal, residuum: %e, multiplier norm: %e",
                    stationarity_residuum,
                    multiplier_norm);

    return false;
  }

  if(slackness_residuum >= slack_eps * (1. + multiplier_norm))
  {
    sleqp_log_debug("Iterate does not satisfy complementary slackness, residuum: %e, multiplier norm: %e",
                    slackness_residuum,
                    multiplier_norm);
    return false;
  }

  return true;
}

SLEQP_RETCODE sleqp_iterate_copy(SleqpIterate* source,
                                 SleqpIterate* target)
{
  SLEQP_CALL(sleqp_sparse_vector_copy(source->primal, target->primal));

  target->func_val = source->func_val;

  SLEQP_CALL(sleqp_sparse_vector_copy(source->func_grad,
                                      target->func_grad));

  SLEQP_CALL(sleqp_sparse_vector_copy(source->cons_val,
                                      target->cons_val));

  SLEQP_CALL(sleqp_sparse_matrix_copy(source->cons_jac,
                                      target->cons_jac));

  SLEQP_CALL(sleqp_working_set_copy(source->working_set,
                                   target->working_set));

  SLEQP_CALL(sleqp_sparse_vector_copy(source->cons_dual,
                                      target->cons_dual));

  SLEQP_CALL(sleqp_sparse_vector_copy(source->vars_dual,
                                      target->vars_dual));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE iterate_free(SleqpIterate** star)
{
  SleqpIterate* iterate = *star;

  if(!iterate)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_sparse_vector_free(&iterate->vars_dual));
  SLEQP_CALL(sleqp_sparse_vector_free(&iterate->cons_dual));

  SLEQP_CALL(sleqp_working_set_release(&iterate->working_set));

  SLEQP_CALL(sleqp_sparse_matrix_release(&iterate->cons_jac));

  SLEQP_CALL(sleqp_sparse_vector_free(&iterate->cons_val));
  SLEQP_CALL(sleqp_sparse_vector_free(&iterate->func_grad));

  SLEQP_CALL(sleqp_sparse_vector_free(&iterate->primal));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_iterate_capture(SleqpIterate* iterate)
{
  ++iterate->refcount;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_iterate_release(SleqpIterate** star)
{
  SleqpIterate* iterate = *star;

  if(!iterate)
  {
    return SLEQP_OKAY;
  }

  if(--iterate->refcount == 0)
  {
    SLEQP_CALL(iterate_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
