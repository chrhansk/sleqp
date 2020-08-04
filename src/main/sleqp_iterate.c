#include "sleqp_iterate.h"

#include <math.h>

#include "sleqp_cmp.h"
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

  assert(sleqp_sparse_vector_valid(x));

  SleqpIterate* iterate = *star;

  *iterate = (SleqpIterate){0};

  iterate->refcount = 1;

  int num_variables = problem->num_variables;
  int num_constraints = problem->num_constraints;

  SLEQP_CALL(sleqp_sparse_vector_create(&iterate->primal,
                                        num_variables,
                                        x->nnz));

  SLEQP_CALL(sleqp_sparse_vector_copy(x, iterate->primal));

  SLEQP_CALL(sleqp_sparse_vector_create(&iterate->func_grad,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&iterate->cons_val,
                                        num_constraints,
                                        0));

  SLEQP_CALL(sleqp_sparse_matrix_create(&iterate->cons_jac,
                                        num_constraints,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_working_set_create(&iterate->working_set,
                                      problem));

  SLEQP_CALL(sleqp_sparse_vector_create(&iterate->cons_dual,
                                        num_constraints,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&iterate->vars_dual,
                                        num_variables,
                                        0));

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

double sleqp_iterate_slackness_residuum(SleqpIterate* iterate,
                                        SleqpProblem* problem)
{
  double residuum = 0.;

  {
    SleqpSparseVec* x = iterate->primal;

    SleqpSparseVec* lb = problem->var_lb;
    SleqpSparseVec* ub = problem->var_ub;

    SleqpSparseVec* d = iterate->vars_dual;

    {
      int k_x = 0, k_ub = 0, k_d = 0;

      int dim = x->dim;

      while((k_x < x->nnz || k_ub < ub->nnz) && k_d < d->nnz)
      {
        bool valid_x = k_x < x->nnz;
        bool valid_ub = k_ub < ub->nnz;
        bool valid_d = k_d < d->nnz;

        int i_x = valid_x ? x->indices[k_x] : dim + 1;
        int i_ub = valid_ub ? ub->indices[k_ub] : dim + 1;
        int i_d = valid_d ? d->indices[k_d] : dim + 1;

        int i_combined = SLEQP_MIN(i_x, i_ub);
        i_combined = SLEQP_MIN(i_combined, i_d);

        valid_x = valid_x && (i_x == i_combined);
        valid_ub = valid_ub && (i_ub == i_combined);
        valid_d = valid_d && (i_d == i_combined);

        double x_val = valid_x ? x->data[k_x] : 0.;
        double ub_val = valid_ub ? ub->data[k_ub] : 0.;
        double d_val = valid_d ? d->data[k_d] : 0.;

        double current_residuum = SLEQP_MAX(ub_val - x_val, 0.) * SLEQP_MAX(d_val, 0.);

        residuum = SLEQP_MAX(residuum, current_residuum);

        if(valid_x)
        {
          ++k_x;
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

    {
      int k_x = 0, k_lb = 0, k_d = 0;

      int dim = x->dim;

      while((k_x < x->nnz || k_lb < lb->nnz) && k_d < d->nnz)
      {
        bool valid_x = k_x < x->nnz;
        bool valid_lb = k_lb < lb->nnz;
        bool valid_d = k_d < d->nnz;

        int i_x = valid_x ? x->indices[k_x] : dim + 1;
        int i_lb = valid_lb ? lb->indices[k_lb] : dim + 1;
        int i_d = valid_d ? d->indices[k_d] : dim + 1;

        int i_combined = SLEQP_MIN(i_x, i_lb);
        i_combined = SLEQP_MIN(i_combined, i_d);

        valid_x = valid_x && (i_x == i_combined);
        valid_lb = valid_lb && (i_lb == i_combined);
        valid_d = valid_d && (i_d == i_combined);

        double x_val = valid_x ? x->data[k_x] : 0.;
        double lb_val = valid_lb ? lb->data[k_lb] : 0.;
        double d_val = valid_d ? d->data[k_d] : 0.;

        double current_residuum = SLEQP_MIN(lb_val - x_val, 0.) * SLEQP_MIN(d_val, 0.);

        residuum = SLEQP_MAX(residuum, current_residuum);

        if(valid_x)
        {
          ++k_x;
        }

        if(valid_lb)
        {
          ++k_lb;
        }

        if(valid_d)
        {
          ++k_d;
        }
      }
    }
  }

  {
    SleqpSparseVec* v = iterate->cons_val;

    SleqpSparseVec* lb = problem->cons_lb;
    SleqpSparseVec* ub = problem->cons_ub;

    SleqpSparseVec* d = iterate->cons_dual;

    {
      int k_v = 0, k_ub = 0, k_d = 0;

      int dim = v->dim;

      while((k_v < v->nnz || k_ub < ub->nnz) && k_d < d->nnz)
      {
        bool valid_v = k_v < v->nnz;
        bool valid_ub = k_ub < ub->nnz;
        bool valid_d = k_d < d->nnz;

        int i_v = valid_v ? v->indices[k_v] : dim + 1;
        int i_ub = valid_ub ? ub->indices[k_ub] : dim + 1;
        int i_d = valid_d ? d->indices[k_d] : dim + 1;

        int i_combined = SLEQP_MIN(i_v, i_ub);
        i_combined = SLEQP_MIN(i_combined, i_d);

        valid_v = valid_v && (i_v == i_combined);
        valid_ub = valid_ub && (i_ub == i_combined);
        valid_d = valid_d && (i_d == i_combined);

        double v_val = valid_v ? v->data[k_v] : 0.;
        double ub_val = valid_ub ? ub->data[k_ub] : 0.;
        double d_val = valid_d ? d->data[k_d] : 0.;

        double current_residuum = SLEQP_MAX(ub_val - v_val, 0.) * SLEQP_MAX(d_val, 0.);

        residuum = SLEQP_MAX(residuum, current_residuum);

        if(valid_v)
        {
          ++k_v;
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


    {
      int k_v = 0, k_lb = 0, k_d = 0;

      int dim = v->dim;

      while((k_v < v->nnz || k_lb < lb->nnz) && k_d < d->nnz)
      {
        bool valid_v = k_v < v->nnz;
        bool valid_lb = k_lb < lb->nnz;
        bool valid_d = k_d < d->nnz;

        int i_v = valid_v ? v->indices[k_v] : dim + 1;
        int i_lb = valid_lb ? lb->indices[k_lb] : dim + 1;
        int i_d = valid_d ? d->indices[k_d] : dim + 1;

        int i_combined = SLEQP_MIN(i_v, i_lb);
        i_combined = SLEQP_MIN(i_combined, i_d);

        valid_v = valid_v && (i_v == i_combined);
        valid_lb = valid_lb && (i_lb == i_combined);
        valid_d = valid_d && (i_d == i_combined);

        double v_val = valid_v ? v->data[k_v] : 0.;
        double lb_val = valid_lb ? lb->data[k_lb] : 0.;
        double d_val = valid_d ? d->data[k_d] : 0.;

        double current_residuum = SLEQP_MIN(lb_val - v_val, 0.) * SLEQP_MIN(d_val, 0.);

        residuum = SLEQP_MAX(residuum, current_residuum);

        if(valid_v)
        {
          ++k_v;
        }

        if(valid_lb)
        {
          ++k_lb;
        }

        if(valid_d)
        {
          ++k_d;
        }
      }
    }
  }

  return residuum;
}

double sleqp_iterate_feasibility_residuum(SleqpIterate* iterate,
                                          SleqpProblem* problem)
{
  SleqpSparseVec* c = iterate->cons_val;
  SleqpSparseVec* lb = problem->cons_lb;
  SleqpSparseVec* ub = problem->cons_ub;

  int k_c = 0, k_lb = 0, k_ub = 0;

  int dim = c->dim;

  double violation = 0.;

  while(k_c < c->nnz || k_lb < lb->nnz || k_ub < ub->nnz)
  {
    bool valid_c = k_c < c->nnz;
    bool valid_lb = k_lb < lb->nnz;
    bool valid_ub = k_ub < ub->nnz;

    int i_c = valid_c ? c->indices[k_c] : dim + 1;
    int i_lb = valid_lb ? lb->indices[k_lb] : dim + 1;
    int i_ub = valid_ub ? ub->indices[k_ub] : dim + 1;

    int i_combined = SLEQP_MIN(i_lb, i_ub);
    i_combined = SLEQP_MIN(i_combined, i_c);

    valid_c = valid_c && (i_c == i_combined);
    valid_lb = valid_lb && (i_lb == i_combined);
    valid_ub = valid_ub && (i_ub == i_combined);

    double c_val = valid_c ? c->data[k_c] : 0.;
    double lb_val = valid_lb ? lb->data[k_lb] : 0.;
    double ub_val = valid_ub ? ub->data[k_ub] : 0.;

    {
      double current_violation = SLEQP_MAX(c_val - ub_val, 0.);

      violation = SLEQP_MAX(violation, current_violation);
    }

    {
      double current_violation = SLEQP_MAX(lb_val - c_val, 0.);

      violation = SLEQP_MAX(violation, current_violation);
    }

    if(valid_c)
    {
      ++k_c;
    }

    if(valid_lb)
    {
      ++k_lb;
    }

    if(valid_ub)
    {
      ++k_ub;
    }
  }

  return violation;
}

SLEQP_RETCODE sleqp_iterate_get_violated_constraints(SleqpIterate* iterate,
                                                     SleqpProblem* problem,
                                                     double tolerance,
                                                     int* violated_constraints,
                                                     int* num_violated_constraints)
{
  double value_norm = 0.;

  {
    value_norm = sleqp_sparse_vector_normsq(iterate->primal);

    value_norm = sqrt(value_norm);
  }

  SleqpSparseVec* c = iterate->cons_val;
  SleqpSparseVec* lb = problem->cons_lb;
  SleqpSparseVec* ub = problem->cons_ub;

  int k_c = 0, k_lb = 0, k_ub = 0;

  int dim = c->dim;

  *num_violated_constraints = 0;

  while(k_c < c->nnz || k_lb < lb->nnz || k_ub < ub->nnz)
  {
    bool valid_c = k_c < c->nnz;
    bool valid_lb = k_lb < lb->nnz;
    bool valid_ub = k_ub < ub->nnz;

    int i_c = valid_c ? c->indices[k_c] : dim + 1;
    int i_lb = valid_lb ? lb->indices[k_lb] : dim + 1;
    int i_ub = valid_ub ? ub->indices[k_ub] : dim + 1;

    int i_combined = SLEQP_MIN(i_lb, i_ub);
    i_combined = SLEQP_MIN(i_combined, i_c);

    valid_c = valid_c && (i_c == i_combined);
    valid_lb = valid_lb && (i_lb == i_combined);
    valid_ub = valid_ub && (i_ub == i_combined);

    double c_val = valid_c ? c->data[k_c] : 0.;
    double lb_val = valid_lb ? lb->data[k_lb] : 0.;
    double ub_val = valid_ub ? ub->data[k_ub] : 0.;

    double lower_violation = SLEQP_MAX(c_val - ub_val, 0.);

    double upper_violation = SLEQP_MAX(lb_val - c_val, 0.);

    double current_violation = SLEQP_MAX(lower_violation, upper_violation);

    if(current_violation > tolerance * (1. + value_norm))
    {
      violated_constraints[(*num_violated_constraints)++] = i_combined;
    }

    if(valid_c)
    {
      ++k_c;
    }

    if(valid_lb)
    {
      ++k_lb;
    }

    if(valid_ub)
    {
      ++k_ub;
    }
  }

  return SLEQP_OKAY;
}

double sleqp_iterate_stationarity_residuum(SleqpIterate* iterate,
                                           SleqpProblem* problem,
                                           double* cache)
{
  double residuum = 0.;

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


  for(int j = 0; j < num_variables; ++j)
  {
    residuum = SLEQP_MAX(residuum, SLEQP_ABS(cache[j]));
  }

  return residuum;
}

bool sleqp_iterate_is_feasible(SleqpIterate* iterate,
                               double feasibility_residuum,
                               double tolerance)
{
  double value_norm = 0.;

  {
    value_norm = sleqp_sparse_vector_normsq(iterate->primal);

    value_norm = sqrt(value_norm);
  }

  if(feasibility_residuum > tolerance * (1. + value_norm))
  {
    sleqp_log_debug("Iterate does not satisfy feasibility, residuum: %e, iterate norm: %e",
                    feasibility_residuum,
                    value_norm);

    return false;
  }

  return true;
}

bool sleqp_iterate_is_optimal(SleqpIterate* iterate,
                              double feasibility_residuum,
                              double slackness_residuum,
                              double stationarity_residuum,
                              double tolerance)
{
  if(!sleqp_iterate_is_feasible(iterate, feasibility_residuum, tolerance))
  {
    return false;
  }

  double multiplier_norm = 0.;

  {
    multiplier_norm += sleqp_sparse_vector_normsq(iterate->cons_dual);
    multiplier_norm += sleqp_sparse_vector_normsq(iterate->vars_dual);

    multiplier_norm = sqrt(multiplier_norm);
  }

  sleqp_log_info("Checking optimality, tolerance: %e", tolerance);

  if(stationarity_residuum >= tolerance * (1. + multiplier_norm))
  {
    sleqp_log_info("Iterate is not optimal, residuum: %e, multiplier norm: %e",
                    stationarity_residuum,
                    multiplier_norm);

    return false;
  }

  if(slackness_residuum >= tolerance * (1. + multiplier_norm))
  {
    sleqp_log_info("Iterate does not satisfy complementary slackness, residuum: %e, multiplier norm: %e",
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
