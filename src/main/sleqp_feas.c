#include "sleqp_feas.h"
#include "sleqp_cmp.h"

SLEQP_RETCODE sleqp_violated_constraint_multipliers(SleqpProblem* problem,
                                                    SleqpSparseVec* cons_vals,
                                                    SleqpSparseVec* multipliers,
                                                    SleqpWorkingSet* working_set,
                                                    double feas_eps)
{
  SleqpSparseVec* lb = problem->cons_lb;
  SleqpSparseVec* ub = problem->cons_ub;
  SleqpSparseVec* v = cons_vals;

  SLEQP_CALL(sleqp_sparse_vector_clear(multipliers));

  SLEQP_CALL(sleqp_sparse_vector_reserve(multipliers, problem->num_constraints));

  const int dim = v->dim;

  assert(lb->dim == dim);
  assert(ub->dim == dim);
  assert(multipliers->dim == dim);

  int k_c = 0, k_lb = 0, k_ub = 0;

  while(k_c < v->nnz || k_lb < lb->nnz || k_ub < ub->nnz)
  {
    double c_val = 0., lb_val = 0., ub_val = 0.;

    bool valid_v = (k_c < v->nnz);
    bool valid_lb = (k_lb < lb->nnz);
    bool valid_ub = (k_ub < ub->nnz);

    int idx = valid_v ? v->indices[k_c] : dim + 1;
    idx = SLEQP_MIN(idx, valid_lb ? lb->indices[k_lb] : dim + 1);
    idx = SLEQP_MIN(idx, valid_ub ? ub->indices[k_ub] : dim + 1);

    if(valid_v && idx == v->indices[k_c])
    {
      c_val = v->data[k_c++];
    }

    if(valid_lb && idx == lb->indices[k_lb])
    {
      lb_val = lb->data[k_lb++];
    }

    if(valid_ub && idx == ub->indices[k_ub])
    {
      ub_val = ub->data[k_ub++];
    }

    if(working_set && sleqp_working_set_get_constraint_state(working_set, idx) != SLEQP_INACTIVE)
    {
      continue;
    }

    const double upper_violation = SLEQP_MAX(c_val - ub_val, 0.);
    const double lower_violation = SLEQP_MAX(lb_val - c_val, 0.);

    if(sleqp_is_pos(upper_violation, feas_eps))
    {
      SLEQP_CALL(sleqp_sparse_vector_push(multipliers,
                                          idx,
                                          1.));
    }
    else if(sleqp_is_pos(lower_violation, feas_eps))
    {
      SLEQP_CALL(sleqp_sparse_vector_push(multipliers,
                                          idx,
                                          -1.));
    }

  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_violated_variable_multipliers(SleqpProblem* problem,
                                                  SleqpSparseVec* primal,
                                                  SleqpSparseVec* multipliers,
                                                  SleqpWorkingSet* working_set,
                                                  double feas_eps)
{
  SleqpSparseVec* lb = problem->var_lb;
  SleqpSparseVec* ub = problem->var_ub;
  SleqpSparseVec* p = primal;

  SLEQP_CALL(sleqp_sparse_vector_clear(multipliers));

  SLEQP_CALL(sleqp_sparse_vector_reserve(multipliers, problem->num_variables));

  const int dim = p->dim;

  assert(lb->dim == dim);
  assert(ub->dim == dim);
  assert(multipliers->dim == dim);

  int k_p = 0, k_lb = 0, k_ub = 0;

  while(k_p < p->nnz || k_lb < lb->nnz || k_ub < ub->nnz)
  {
    double p_val = 0., lb_val = 0., ub_val = 0.;

    bool valid_p = (k_p < p->nnz);
    bool valid_lb = (k_lb < lb->nnz);
    bool valid_ub = (k_ub < ub->nnz);

    int idx = valid_p ? p->indices[k_p] : dim + 1;
    idx = SLEQP_MIN(idx, valid_lb ? lb->indices[k_lb] : dim + 1);
    idx = SLEQP_MIN(idx, valid_ub ? ub->indices[k_ub] : dim + 1);

    if(valid_p && idx == p->indices[k_p])
    {
      p_val = p->data[k_p++];
    }

    if(valid_lb && idx == lb->indices[k_lb])
    {
      lb_val = lb->data[k_lb++];
    }

    if(valid_ub && idx == ub->indices[k_ub])
    {
      ub_val = ub->data[k_ub++];
    }

    if(working_set && sleqp_working_set_get_variable_state(working_set, idx) != SLEQP_INACTIVE)
    {
      continue;
    }

    const double upper_violation = SLEQP_MAX(p_val - ub_val, 0.);
    const double lower_violation = SLEQP_MAX(lb_val - p_val, 0.);

    if(sleqp_is_pos(upper_violation, feas_eps))
    {
      SLEQP_CALL(sleqp_sparse_vector_push(multipliers,
                                          idx,
                                          1.));
    }
    else if(sleqp_is_pos(lower_violation, feas_eps))
    {
      SLEQP_CALL(sleqp_sparse_vector_push(multipliers,
                                          idx,
                                          -1.));
    }

  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_violated_constraints(SleqpProblem* problem,
                                         SleqpSparseVec* cons_val,
                                         int* violated_constraints,
                                         int* num_violated_constraints,
                                         double feas_eps)
{
  SleqpSparseVec* c = cons_val;
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

    const double upper_violation = SLEQP_MAX(c_val - ub_val, 0.);
    const double lower_violation = SLEQP_MAX(lb_val - c_val, 0.);

    if(sleqp_is_pos(upper_violation, feas_eps))
    {
      violated_constraints[(*num_violated_constraints)++] = i_combined;
    }
    else if(sleqp_is_pos(lower_violation, feas_eps))
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

SLEQP_RETCODE sleqp_violation_values(SleqpProblem* problem,
                                     SleqpIterate* iterate,
                                     double zero_eps,
                                     SleqpSparseVec* violation)
{
  int num_constraints = problem->num_constraints;

  SleqpSparseVec* lb = problem->cons_lb;
  SleqpSparseVec* ub = problem->cons_ub;
  SleqpSparseVec* c = sleqp_iterate_get_cons_val(iterate);

  const int dim = c->dim;

  assert(dim == num_constraints);
  assert(dim == lb->dim);
  assert(dim == ub->dim);

  {
    int max_size = lb->nnz + ub->nnz + c->nnz;
    max_size = SLEQP_MIN(max_size, num_constraints);

    SLEQP_CALL(sleqp_sparse_vector_clear(violation));
    SLEQP_CALL(sleqp_sparse_vector_resize(violation, dim));
    SLEQP_CALL(sleqp_sparse_vector_reserve(violation, max_size));
  }

  int k_c = 0, k_lb = 0, k_ub = 0;

  while(k_c < c->nnz || k_lb < lb->nnz || k_ub < ub->nnz)
  {
    double c_val = 0., lb_val = 0., ub_val = 0.;

    bool valid_c = (k_c < c->nnz);
    bool valid_lb = (k_lb < lb->nnz);
    bool valid_ub = (k_ub < ub->nnz);

    int idx = valid_c ? c->indices[k_c] : dim + 1;
    idx = SLEQP_MIN(idx, valid_lb ? lb->indices[k_lb] : dim + 1);
    idx = SLEQP_MIN(idx, valid_ub ? ub->indices[k_ub] : dim + 1);

    if(valid_c && idx == c->indices[k_c])
    {
      c_val = c->data[k_c++];
    }

    if(valid_lb && idx == lb->indices[k_lb])
    {
      lb_val = lb->data[k_lb++];
    }

    if(valid_ub && idx == ub->indices[k_ub])
    {
      ub_val = ub->data[k_ub++];
    }

    const double upper_violation = SLEQP_MAX(c_val - ub_val, 0.);
    const double lower_violation = SLEQP_MAX(lb_val - c_val, 0.);

    if(sleqp_is_pos(upper_violation, zero_eps))
    {
      SLEQP_CALL(sleqp_sparse_vector_push(violation,
                                          idx,
                                          c_val - ub_val));
    }
    else if(sleqp_is_pos(lower_violation, zero_eps))
    {
      SLEQP_CALL(sleqp_sparse_vector_push(violation,
                                          idx,
                                          c_val - lb_val));
    }

  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_violation_inf_norm(SleqpProblem* problem,
                                       SleqpSparseVec* cons_val,
                                       double zero_eps,
                                       double* max_violation)
{
  int num_constraints = problem->num_constraints;

  SleqpSparseVec* lb = problem->cons_lb;
  SleqpSparseVec* ub = problem->cons_ub;
  SleqpSparseVec* c = cons_val;

  const int dim = c->dim;

  assert(dim == num_constraints);
  assert(dim == lb->dim);
  assert(dim == ub->dim);

  (*max_violation) = 0.;

  int k_c = 0, k_lb = 0, k_ub = 0;

  while(k_c < c->nnz || k_lb < lb->nnz || k_ub < ub->nnz)
  {
    double c_val = 0., lb_val = 0., ub_val = 0.;

    bool valid_c = (k_c < c->nnz);
    bool valid_lb = (k_lb < lb->nnz);
    bool valid_ub = (k_ub < ub->nnz);

    int idx = valid_c ? c->indices[k_c] : dim + 1;
    idx = SLEQP_MIN(idx, valid_lb ? lb->indices[k_lb] : dim + 1);
    idx = SLEQP_MIN(idx, valid_ub ? ub->indices[k_ub] : dim + 1);

    if(valid_c && idx == c->indices[k_c])
    {
      c_val = c->data[k_c++];
    }

    if(valid_lb && idx == lb->indices[k_lb])
    {
      lb_val = lb->data[k_lb++];
    }

    if(valid_ub && idx == ub->indices[k_ub])
    {
      ub_val = ub->data[k_ub++];
    }

    const double upper_violation = SLEQP_MAX(c_val - ub_val, 0.);
    const double lower_violation = SLEQP_MAX(lb_val - c_val, 0.);

    if(sleqp_is_pos(upper_violation, zero_eps))
    {
      (*max_violation) = SLEQP_MAX(upper_violation, (*max_violation));
    }
    else if(sleqp_is_pos(lower_violation, zero_eps))
    {
      (*max_violation) = SLEQP_MAX(lower_violation, (*max_violation));
    }

  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_violation_one_norm(SleqpProblem* problem,
                                       SleqpSparseVec* cons_val,
                                       double zero_eps,
                                       double* total_violation)
{
  int num_constraints = problem->num_constraints;

  SleqpSparseVec* lb = problem->cons_lb;
  SleqpSparseVec* ub = problem->cons_ub;
  SleqpSparseVec* c = cons_val;

  const int dim = c->dim;

  assert(dim == num_constraints);
  assert(dim == lb->dim);
  assert(dim == ub->dim);

  (*total_violation) = 0.;

  int k_c = 0, k_lb = 0, k_ub = 0;

  while(k_c < c->nnz || k_lb < lb->nnz || k_ub < ub->nnz)
  {
    double c_val = 0., lb_val = 0., ub_val = 0.;

    bool valid_c = (k_c < c->nnz);
    bool valid_lb = (k_lb < lb->nnz);
    bool valid_ub = (k_ub < ub->nnz);

    int idx = valid_c ? c->indices[k_c] : dim + 1;
    idx = SLEQP_MIN(idx, valid_lb ? lb->indices[k_lb] : dim + 1);
    idx = SLEQP_MIN(idx, valid_ub ? ub->indices[k_ub] : dim + 1);

    if(valid_c && idx == c->indices[k_c])
    {
      c_val = c->data[k_c++];
    }

    if(valid_lb && idx == lb->indices[k_lb])
    {
      lb_val = lb->data[k_lb++];
    }

    if(valid_ub && idx == ub->indices[k_ub])
    {
      ub_val = ub->data[k_ub++];
    }

    const double upper_violation = SLEQP_MAX(c_val - ub_val, 0.);
    const double lower_violation = SLEQP_MAX(lb_val - c_val, 0.);

    if(sleqp_is_pos(upper_violation, zero_eps))
    {
      (*total_violation) += upper_violation;
    }
    else if(sleqp_is_pos(lower_violation, zero_eps))
    {
      (*total_violation) += lower_violation;
    }
  }

  return SLEQP_OKAY;
}
