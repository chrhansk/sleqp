#include "iterate.h"

#include <math.h>

#include "cmp.h"
#include "fail.h"
#include "feas.h"
#include "mem.h"
#include "pub_types.h"
#include "working_set.h"

struct SleqpIterate
{
  int refcount;
  /**
   * The current point. Has dimension = num_variables.
   **/
  SleqpVec* primal;

  /**
   * The current function value
   **/
  double obj_val;

  /**
   * The current function gradient. Has dimension = num_variables.
   **/
  SleqpVec* obj_grad;

  /**
   * The current constraint values. Has dimension = num_constraints.
   **/
  SleqpVec* cons_val;

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
  SleqpVec* cons_dual;

  /**
   * The dual values of the variable bounds. Has dimension = num_variables.
   */
  SleqpVec* vars_dual;
};

SLEQP_RETCODE
sleqp_iterate_create(SleqpIterate** star,
                     SleqpProblem* problem,
                     const SleqpVec* x)
{
  SLEQP_CALL(sleqp_malloc(star));

  sleqp_assert_msg(sleqp_vec_is_valid(x), "Invalid primal values");

  sleqp_assert_msg(sleqp_vec_is_finite(x), "Infinite primal values");

  SleqpIterate* iterate = *star;

  *iterate = (SleqpIterate){0};

  iterate->refcount = 1;

  const int num_variables   = sleqp_problem_num_vars(problem);
  const int num_constraints = sleqp_problem_num_cons(problem);

  SLEQP_CALL(sleqp_vec_create(&iterate->primal, num_variables, x->nnz));

  SLEQP_CALL(sleqp_vec_copy(x, iterate->primal));

  SLEQP_CALL(sleqp_vec_create_empty(&iterate->obj_grad, num_variables));

  SLEQP_CALL(sleqp_vec_create_empty(&iterate->cons_val, num_constraints));

  SLEQP_CALL(sleqp_sparse_matrix_create(&iterate->cons_jac,
                                        num_constraints,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_working_set_create(&iterate->working_set, problem));

  SLEQP_CALL(sleqp_vec_create_empty(&iterate->cons_dual, num_constraints));

  SLEQP_CALL(sleqp_vec_create_empty(&iterate->vars_dual, num_variables));

  return SLEQP_OKAY;
}

SleqpVec*
sleqp_iterate_primal(const SleqpIterate* iterate)
{
  return iterate->primal;
}

double
sleqp_iterate_obj_val(const SleqpIterate* iterate)
{
  return iterate->obj_val;
}

SLEQP_RETCODE
sleqp_iterate_set_obj_val(SleqpIterate* iterate, double value)
{
  iterate->obj_val = value;
  return SLEQP_OKAY;
}

SleqpVec*
sleqp_iterate_obj_grad(const SleqpIterate* iterate)
{
  return iterate->obj_grad;
}

SleqpVec*
sleqp_iterate_cons_val(const SleqpIterate* iterate)
{
  return iterate->cons_val;
}

SleqpSparseMatrix*
sleqp_iterate_cons_jac(const SleqpIterate* iterate)
{
  return iterate->cons_jac;
}

SleqpWorkingSet*
sleqp_iterate_working_set(const SleqpIterate* iterate)
{
  return iterate->working_set;
}

SleqpVec*
sleqp_iterate_cons_dual(const SleqpIterate* iterate)
{
  return iterate->cons_dual;
}

SleqpVec*
sleqp_iterate_vars_dual(const SleqpIterate* iterate)
{
  return iterate->vars_dual;
}

static SLEQP_RETCODE
slack_residuum(const SleqpVec* v,
               const SleqpVec* lb,
               const SleqpVec* ub,
               const SleqpVec* d,
               double* residuum)
{
  const int dim = v->dim;

  *residuum = 0.;

  assert(lb->dim == dim);
  assert(ub->dim == dim);
  assert(d->dim == dim);

  {
    {
      int k_v = 0, k_ub = 0, k_lb = 0, k_d = 0;

      int dim = v->dim;

      while ((k_v < v->nnz || k_lb < lb->nnz || k_ub < ub->nnz) && k_d < d->nnz)
      {
        bool valid_v  = k_v < v->nnz;
        bool valid_lb = k_lb < lb->nnz;
        bool valid_ub = k_ub < ub->nnz;
        bool valid_d  = k_d < d->nnz;

        const int i_v  = valid_v ? v->indices[k_v] : dim + 1;
        const int i_lb = valid_lb ? lb->indices[k_lb] : dim + 1;
        const int i_ub = valid_ub ? ub->indices[k_ub] : dim + 1;
        const int i_d  = valid_d ? d->indices[k_d] : dim + 1;

        int i_combined = SLEQP_MIN(i_lb, i_ub);
        i_combined     = SLEQP_MIN(i_combined, i_d);
        i_combined     = SLEQP_MIN(i_combined, i_v);

        valid_v  = valid_v && (i_v == i_combined);
        valid_lb = valid_lb && (i_lb == i_combined);
        valid_ub = valid_ub && (i_ub == i_combined);
        valid_d  = valid_d && (i_d == i_combined);

        const double v_val  = valid_v ? v->data[k_v] : 0.;
        const double lb_val = valid_lb ? lb->data[k_lb] : 0.;
        const double ub_val = valid_ub ? ub->data[k_ub] : 0.;
        const double d_val  = valid_d ? d->data[k_d] : 0.;

        double current_residuum = 0.;

        // use signs of dual variables with respect to working set
        if (d_val >= 0.)
        {
          current_residuum = SLEQP_MAX(ub_val - v_val, 0.) * d_val;
        }
        else
        {
          current_residuum = SLEQP_MAX(v_val - lb_val, 0.) * d_val;
        }

        *residuum = SLEQP_MAX(*residuum, SLEQP_ABS(current_residuum));

        if (valid_v)
        {
          ++k_v;
        }

        if (valid_lb)
        {
          ++k_lb;
        }

        if (valid_ub)
        {
          ++k_ub;
        }

        if (valid_d)
        {
          ++k_d;
        }
      }
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_iterate_slackness_residuum(SleqpProblem* problem,
                                 SleqpIterate* iterate,
                                 double* slackness_residuum)
{
  (*slackness_residuum) = 0.;

  double current_residuum;

  SLEQP_CALL(slack_residuum(sleqp_iterate_cons_val(iterate),
                            sleqp_problem_cons_lb(problem),
                            sleqp_problem_cons_ub(problem),
                            sleqp_iterate_cons_dual(iterate),
                            &current_residuum));

  (*slackness_residuum) = SLEQP_MAX((*slackness_residuum), current_residuum);

  SLEQP_CALL(slack_residuum(sleqp_iterate_primal(iterate),
                            sleqp_problem_vars_lb(problem),
                            sleqp_problem_vars_ub(problem),
                            sleqp_iterate_vars_dual(iterate),
                            &current_residuum));

  (*slackness_residuum) = SLEQP_MAX((*slackness_residuum), current_residuum);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
slack_residuals(const SleqpVec* v,
                const SleqpVec* lb,
                const SleqpVec* ub,
                const SleqpVec* d,
                SleqpVec* r,
                double zero_eps)
{
  const int dim = v->dim;

  SLEQP_CALL(sleqp_vec_clear(r));

  SLEQP_CALL(sleqp_vec_resize(r, dim));
  SLEQP_CALL(sleqp_vec_reserve(r, dim));

  {
    {
      int k_v = 0, k_ub = 0, k_lb = 0, k_d = 0;

      int dim = v->dim;

      while ((k_v < v->nnz || k_lb < lb->nnz || k_ub < ub->nnz) && k_d < d->nnz)
      {
        bool valid_v  = k_v < v->nnz;
        bool valid_lb = k_lb < lb->nnz;
        bool valid_ub = k_ub < ub->nnz;
        bool valid_d  = k_d < d->nnz;

        const int i_v  = valid_v ? v->indices[k_v] : dim + 1;
        const int i_lb = valid_lb ? lb->indices[k_lb] : dim + 1;
        const int i_ub = valid_ub ? ub->indices[k_ub] : dim + 1;
        const int i_d  = valid_d ? d->indices[k_d] : dim + 1;

        int i_combined = SLEQP_MIN(i_lb, i_ub);
        i_combined     = SLEQP_MIN(i_combined, i_d);
        i_combined     = SLEQP_MIN(i_combined, i_v);

        valid_v  = valid_v && (i_v == i_combined);
        valid_lb = valid_lb && (i_lb == i_combined);
        valid_ub = valid_ub && (i_ub == i_combined);
        valid_d  = valid_d && (i_d == i_combined);

        const double v_val  = valid_v ? v->data[k_v] : 0.;
        const double lb_val = valid_lb ? lb->data[k_lb] : 0.;
        const double ub_val = valid_ub ? ub->data[k_ub] : 0.;
        const double d_val  = valid_d ? d->data[k_d] : 0.;

        double current_residuum = 0.;

        // use signs of dual variables with respect to working set
        if (d_val >= 0.)
        {
          current_residuum = SLEQP_MAX(ub_val - v_val, 0.) * d_val;
        }
        else
        {
          current_residuum = SLEQP_MAX(v_val - lb_val, 0.) * d_val;
        }

        if (!sleqp_is_zero(current_residuum, zero_eps))
        {
          SLEQP_CALL(sleqp_vec_push(r, i_combined, current_residuum));
        }

        if (valid_v)
        {
          ++k_v;
        }

        if (valid_lb)
        {
          ++k_lb;
        }

        if (valid_ub)
        {
          ++k_ub;
        }

        if (valid_d)
        {
          ++k_d;
        }
      }
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_iterate_vars_slackness_residuals(SleqpProblem* problem,
                                       SleqpIterate* iterate,
                                       SleqpVec* residuals,
                                       double zero_eps)
{
  SLEQP_CALL(slack_residuals(sleqp_iterate_primal(iterate),
                             sleqp_problem_vars_lb(problem),
                             sleqp_problem_vars_ub(problem),
                             sleqp_iterate_vars_dual(iterate),
                             residuals,
                             zero_eps));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_iterate_cons_slackness_residuals(SleqpProblem* problem,
                                       SleqpIterate* iterate,
                                       SleqpVec* residuals,
                                       double zero_eps)
{
  SLEQP_CALL(slack_residuals(sleqp_iterate_cons_val(iterate),
                             sleqp_problem_cons_lb(problem),
                             sleqp_problem_cons_ub(problem),
                             sleqp_iterate_cons_dual(iterate),
                             residuals,
                             zero_eps));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_iterate_feasibility_residuum(SleqpProblem* problem,
                                   SleqpIterate* iterate,
                                   double* feasibility_residuum)
{
  SLEQP_CALL(
    sleqp_max_violation(problem, iterate->cons_val, feasibility_residuum));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_iterate_get_violated_constraints(SleqpProblem* problem,
                                       SleqpIterate* iterate,
                                       int* violated_constraints,
                                       int* num_violated_constraints,
                                       double feas_eps)
{
  SLEQP_CALL(sleqp_violated_constraints(problem,
                                        iterate->cons_val,
                                        violated_constraints,
                                        num_violated_constraints));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
write_stationarity_resiudals_to_cache(SleqpIterate* iterate,
                                      SleqpProblem* problem,
                                      double* cache)
{
  const int num_variables   = sleqp_problem_num_vars(problem);
  const int num_constraints = sleqp_problem_num_cons(problem);

  const SleqpSparseMatrix* cons_jac = sleqp_iterate_cons_jac(iterate);
  const SleqpVec* obj_grad          = sleqp_iterate_obj_grad(iterate);

  const SleqpVec* cons_dual = sleqp_iterate_cons_dual(iterate);
  const SleqpVec* vars_dual = sleqp_iterate_vars_dual(iterate);

  const int num_rows = sleqp_sparse_matrix_num_rows(cons_jac);
  const int num_cols = sleqp_sparse_matrix_num_cols(cons_jac);

  const int* cons_jac_cols    = sleqp_sparse_matrix_cols(cons_jac);
  const int* cons_jac_rows    = sleqp_sparse_matrix_rows(cons_jac);
  const double* cons_jac_data = sleqp_sparse_matrix_data(cons_jac);

  assert(num_variables == num_cols);
  assert(num_constraints == num_rows);

  for (int j = 0; j < num_variables; ++j)
  {
    int k_d = 0;
    int k_j = cons_jac_cols[j];

    double sum = 0.;

    while (k_d < cons_dual->nnz && k_j < cons_jac_cols[j + 1])
    {
      const int d_idx = cons_dual->indices[k_d];
      const int j_idx = cons_jac_rows[k_j];

      if (d_idx < j_idx)
      {
        ++k_d;
      }
      else if (d_idx > j_idx)
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

  for (int k = 0; k < vars_dual->nnz; ++k)
  {
    cache[vars_dual->indices[k]] += vars_dual->data[k];
  }

  for (int k = 0; k < obj_grad->nnz; ++k)
  {
    cache[obj_grad->indices[k]] += obj_grad->data[k];
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_iterate_stationarity_residuals(SleqpProblem* problem,
                                     SleqpIterate* iterate,
                                     double* cache,
                                     SleqpVec* residuals,
                                     double zero_eps)
{
  const int num_variables = sleqp_problem_num_vars(problem);

  SLEQP_CALL(write_stationarity_resiudals_to_cache(iterate, problem, cache));

  SLEQP_CALL(sleqp_vec_set_from_raw(residuals, cache, num_variables, zero_eps));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_iterate_stationarity_residuum(SleqpProblem* problem,
                                    SleqpIterate* iterate,
                                    double* cache,
                                    double* stationarity_residuum)
{
  const int num_variables = sleqp_problem_num_vars(problem);

  SLEQP_CALL(write_stationarity_resiudals_to_cache(iterate, problem, cache));

  (*stationarity_residuum) = 0.;

  for (int j = 0; j < num_variables; ++j)
  {
    (*stationarity_residuum)
      = SLEQP_MAX((*stationarity_residuum), SLEQP_ABS(cache[j]));
  }

  return SLEQP_OKAY;
}

bool
sleqp_iterate_is_feasible(SleqpIterate* iterate,
                          double feasibility_residuum,
                          double feasibility_tolerance)
{
  return feasibility_residuum <= feasibility_tolerance;
}

bool
sleqp_iterate_is_optimal(SleqpIterate* iterate,
                         SleqpParams* params,
                         double feasibility_residuum,
                         double slackness_residuum,
                         double stationarity_residuum)
{
  const double feas_eps = sleqp_params_value(params, SLEQP_PARAM_FEAS_TOL);

  const double slack_eps = sleqp_params_value(params, SLEQP_PARAM_SLACK_TOL);

  const double stat_eps = sleqp_params_value(params, SLEQP_PARAM_STAT_TOL);

  if (!sleqp_iterate_is_feasible(iterate, feasibility_residuum, feas_eps))
  {
    sleqp_log_debug("Iterate is not feasible, residuum: %e",
                    feasibility_residuum);

    return false;
  }

  if (stationarity_residuum >= stat_eps)
  {
    sleqp_log_debug("Iterate is not optimal, stationarity residuum: %e",
                    stationarity_residuum);

    return false;
  }

  if (slackness_residuum >= slack_eps)
  {
    sleqp_log_debug("Iterate is not optimal, slackness residuum: %e",
                    slackness_residuum);
    return false;
  }

  return true;
}

SLEQP_RETCODE
sleqp_iterate_reserve(SleqpIterate* iterate, SleqpProblem* problem)
{
  int obj_grad_nnz  = SLEQP_NONE;
  int cons_val_nnz  = SLEQP_NONE;
  int cons_jac_nnz  = SLEQP_NONE;
  int hess_prod_nnz = SLEQP_NONE;

  SLEQP_CALL(sleqp_problem_nonzeros(problem,
                                    &obj_grad_nnz,
                                    &cons_val_nnz,
                                    &cons_jac_nnz,
                                    &hess_prod_nnz));

  const int num_vars = sleqp_problem_num_vars(problem);
  const int num_cons = sleqp_problem_num_cons(problem);

  if (obj_grad_nnz == SLEQP_NONE)
  {
    obj_grad_nnz = num_vars;
  }

  if (cons_val_nnz == SLEQP_NONE)
  {
    cons_val_nnz = num_cons;
  }

  if (cons_jac_nnz == SLEQP_NONE)
  {
    cons_jac_nnz = num_vars * num_cons;
  }

  SLEQP_CALL(sleqp_vec_reserve(sleqp_iterate_obj_grad(iterate), obj_grad_nnz));
  SLEQP_CALL(sleqp_vec_reserve(sleqp_iterate_cons_val(iterate), cons_val_nnz));

  SLEQP_CALL(
    sleqp_sparse_matrix_reserve(sleqp_iterate_cons_jac(iterate), cons_jac_nnz));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_iterate_copy(const SleqpIterate* source, SleqpIterate* target)
{
  SLEQP_CALL(sleqp_vec_copy(source->primal, target->primal));

  target->obj_val = source->obj_val;

  SLEQP_CALL(sleqp_vec_copy(source->obj_grad, target->obj_grad));

  SLEQP_CALL(sleqp_vec_copy(source->cons_val, target->cons_val));

  SLEQP_CALL(sleqp_sparse_matrix_copy(source->cons_jac, target->cons_jac));

  SLEQP_CALL(sleqp_working_set_copy(source->working_set, target->working_set));

  SLEQP_CALL(sleqp_vec_copy(source->cons_dual, target->cons_dual));

  SLEQP_CALL(sleqp_vec_copy(source->vars_dual, target->vars_dual));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
iterate_free(SleqpIterate** star)
{
  SleqpIterate* iterate = *star;

  if (!iterate)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_vec_free(&iterate->vars_dual));
  SLEQP_CALL(sleqp_vec_free(&iterate->cons_dual));

  SLEQP_CALL(sleqp_working_set_release(&iterate->working_set));

  SLEQP_CALL(sleqp_sparse_matrix_release(&iterate->cons_jac));

  SLEQP_CALL(sleqp_vec_free(&iterate->cons_val));
  SLEQP_CALL(sleqp_vec_free(&iterate->obj_grad));

  SLEQP_CALL(sleqp_vec_free(&iterate->primal));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_iterate_capture(SleqpIterate* iterate)
{
  ++iterate->refcount;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_iterate_release(SleqpIterate** star)
{
  SleqpIterate* iterate = *star;

  if (!iterate)
  {
    return SLEQP_OKAY;
  }

  if (--iterate->refcount == 0)
  {
    SLEQP_CALL(iterate_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
