#include "feas.h"

#include "cmp.h"
#include "log.h"

SLEQP_RETCODE
sleqp_violated_constraint_multipliers(SleqpProblem* problem,
                                      const SleqpVec* cons_vals,
                                      SleqpVec* multipliers,
                                      SleqpWorkingSet* working_set)
{
  const SleqpVec* lb = sleqp_problem_cons_lb(problem);
  const SleqpVec* ub = sleqp_problem_cons_ub(problem);
  const SleqpVec* v  = cons_vals;

  const int num_constraints = sleqp_problem_num_cons(problem);

  SLEQP_CALL(sleqp_vec_clear(multipliers));

  SLEQP_CALL(sleqp_vec_reserve(multipliers, num_constraints));

  const int dim = v->dim;

  assert(lb->dim == dim);
  assert(ub->dim == dim);
  assert(multipliers->dim == dim);

  int k_c = 0, k_lb = 0, k_ub = 0;

  while (k_c < v->nnz || k_lb < lb->nnz || k_ub < ub->nnz)
  {
    double c_val = 0., lb_val = 0., ub_val = 0.;

    bool valid_v  = (k_c < v->nnz);
    bool valid_lb = (k_lb < lb->nnz);
    bool valid_ub = (k_ub < ub->nnz);

    int idx = valid_v ? v->indices[k_c] : dim + 1;
    idx     = SLEQP_MIN(idx, valid_lb ? lb->indices[k_lb] : dim + 1);
    idx     = SLEQP_MIN(idx, valid_ub ? ub->indices[k_ub] : dim + 1);

    if (valid_v && idx == v->indices[k_c])
    {
      c_val = v->data[k_c++];
    }

    if (valid_lb && idx == lb->indices[k_lb])
    {
      lb_val = lb->data[k_lb++];
    }

    if (valid_ub && idx == ub->indices[k_ub])
    {
      ub_val = ub->data[k_ub++];
    }

    if (working_set
        && sleqp_working_set_cons_state(working_set, idx) != SLEQP_INACTIVE)
    {
      continue;
    }

    const double upper_violation = SLEQP_MAX(c_val - ub_val, 0.);
    const double lower_violation = SLEQP_MAX(lb_val - c_val, 0.);

    if (upper_violation > 0.)
    {
      SLEQP_CALL(sleqp_vec_push(multipliers, idx, 1.));
    }
    else if (lower_violation > 0.)
    {
      SLEQP_CALL(sleqp_vec_push(multipliers, idx, -1.));
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_violated_constraints(SleqpProblem* problem,
                           SleqpVec* cons_val,
                           int* violated_constraints,
                           int* num_violated_constraints)
{
  SleqpVec* c        = cons_val;
  const SleqpVec* lb = sleqp_problem_cons_lb(problem);
  const SleqpVec* ub = sleqp_problem_cons_ub(problem);

  int k_c = 0, k_lb = 0, k_ub = 0;

  int dim = c->dim;

  *num_violated_constraints = 0;

  while (k_c < c->nnz || k_lb < lb->nnz || k_ub < ub->nnz)
  {
    bool valid_c  = k_c < c->nnz;
    bool valid_lb = k_lb < lb->nnz;
    bool valid_ub = k_ub < ub->nnz;

    int i_c  = valid_c ? c->indices[k_c] : dim + 1;
    int i_lb = valid_lb ? lb->indices[k_lb] : dim + 1;
    int i_ub = valid_ub ? ub->indices[k_ub] : dim + 1;

    int i_combined = SLEQP_MIN(i_lb, i_ub);
    i_combined     = SLEQP_MIN(i_combined, i_c);

    valid_c  = valid_c && (i_c == i_combined);
    valid_lb = valid_lb && (i_lb == i_combined);
    valid_ub = valid_ub && (i_ub == i_combined);

    double c_val  = valid_c ? c->data[k_c] : 0.;
    double lb_val = valid_lb ? lb->data[k_lb] : 0.;
    double ub_val = valid_ub ? ub->data[k_ub] : 0.;

    const double upper_violation = SLEQP_MAX(c_val - ub_val, 0.);
    const double lower_violation = SLEQP_MAX(lb_val - c_val, 0.);

    if (upper_violation > 0.)
    {
      violated_constraints[(*num_violated_constraints)++] = i_combined;
    }
    else if (lower_violation > 0.)
    {
      violated_constraints[(*num_violated_constraints)++] = i_combined;
    }

    if (valid_c)
    {
      ++k_c;
    }

    if (valid_lb)
    {
      ++k_lb;
    }

    if (valid_ub)
    {
      ++k_ub;
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_violation_values(SleqpProblem* problem,
                       const SleqpVec* cons_val,
                       SleqpVec* violation)
{
  const int num_constraints = sleqp_problem_num_cons(problem);

  const SleqpVec* lb = sleqp_problem_cons_lb(problem);
  const SleqpVec* ub = sleqp_problem_cons_ub(problem);
  const SleqpVec* c  = cons_val;

  const int dim = c->dim;

  assert(dim == num_constraints);
  assert(dim == lb->dim);
  assert(dim == ub->dim);

  {
    int max_size = lb->nnz + ub->nnz + c->nnz;
    max_size     = SLEQP_MIN(max_size, num_constraints);

    SLEQP_CALL(sleqp_vec_clear(violation));
    SLEQP_CALL(sleqp_vec_resize(violation, dim));
    SLEQP_CALL(sleqp_vec_reserve(violation, max_size));
  }

  int k_c = 0, k_lb = 0, k_ub = 0;

  while (k_c < c->nnz || k_lb < lb->nnz || k_ub < ub->nnz)
  {
    double c_val = 0., lb_val = 0., ub_val = 0.;

    bool valid_c  = (k_c < c->nnz);
    bool valid_lb = (k_lb < lb->nnz);
    bool valid_ub = (k_ub < ub->nnz);

    int idx = valid_c ? c->indices[k_c] : dim + 1;
    idx     = SLEQP_MIN(idx, valid_lb ? lb->indices[k_lb] : dim + 1);
    idx     = SLEQP_MIN(idx, valid_ub ? ub->indices[k_ub] : dim + 1);

    if (valid_c && idx == c->indices[k_c])
    {
      c_val = c->data[k_c++];
    }

    if (valid_lb && idx == lb->indices[k_lb])
    {
      lb_val = lb->data[k_lb++];
    }

    if (valid_ub && idx == ub->indices[k_ub])
    {
      ub_val = ub->data[k_ub++];
    }

    const double upper_violation = SLEQP_MAX(c_val - ub_val, 0.);
    const double lower_violation = SLEQP_MAX(lb_val - c_val, 0.);

    if (upper_violation > 0.)
    {
      SLEQP_CALL(sleqp_vec_push(violation, idx, c_val - ub_val));
    }
    else if (lower_violation > 0.)
    {
      SLEQP_CALL(sleqp_vec_push(violation, idx, c_val - lb_val));
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
feasibility_residuals(SleqpProblem* problem,
                      const SleqpVec* cons_val,
                      SleqpVec* residuals,
                      SleqpWorkingSet* working_set,
                      bool signed_residuals)
{
  const int num_constraints = sleqp_problem_num_cons(problem);

  const SleqpVec* lb = sleqp_problem_cons_lb(problem);
  const SleqpVec* ub = sleqp_problem_cons_ub(problem);
  const SleqpVec* c  = cons_val;

  const int dim = c->dim;

  assert(dim == num_constraints);
  assert(dim == lb->dim);
  assert(dim == ub->dim);

  {
    int max_size = lb->nnz + ub->nnz + c->nnz;
    max_size     = SLEQP_MIN(max_size, num_constraints);

    SLEQP_CALL(sleqp_vec_clear(residuals));
    SLEQP_CALL(sleqp_vec_resize(residuals, dim));
    SLEQP_CALL(sleqp_vec_reserve(residuals, max_size));
  }

  int k_c = 0, k_lb = 0, k_ub = 0;

  while (k_c < c->nnz || k_lb < lb->nnz || k_ub < ub->nnz)
  {
    double c_val = 0., lb_val = 0., ub_val = 0.;

    bool valid_c  = (k_c < c->nnz);
    bool valid_lb = (k_lb < lb->nnz);
    bool valid_ub = (k_ub < ub->nnz);

    int idx = valid_c ? c->indices[k_c] : dim + 1;
    idx     = SLEQP_MIN(idx, valid_lb ? lb->indices[k_lb] : dim + 1);
    idx     = SLEQP_MIN(idx, valid_ub ? ub->indices[k_ub] : dim + 1);

    if (valid_c && idx == c->indices[k_c])
    {
      c_val = c->data[k_c++];
    }

    if (valid_lb && idx == lb->indices[k_lb])
    {
      lb_val = lb->data[k_lb++];
    }

    if (valid_ub && idx == ub->indices[k_ub])
    {
      ub_val = ub->data[k_ub++];
    }

    if (working_set
        && sleqp_working_set_cons_state(working_set, idx) != SLEQP_INACTIVE)
    {
      continue;
    }

    const double upper_residual = SLEQP_MAX(c_val - ub_val, 0.);
    const double lower_residual = SLEQP_MAX(lb_val - c_val, 0.);

    assert(upper_residual == 0. || lower_residual == 0.);

    double residual = 0.;

    if (signed_residuals)
    {
      residual = (upper_residual != 0.) ? upper_residual : (-lower_residual);
    }
    else
    {
      residual = SLEQP_MAX(lower_residual, upper_residual);
    }

    if (residual != 0.)
    {
      SLEQP_CALL(sleqp_vec_push(residuals, idx, residual));
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_feasibility_residuals(SleqpProblem* problem,
                            const SleqpVec* cons_val,
                            SleqpVec* residuals,
                            SleqpWorkingSet* working_set)
{
  return feasibility_residuals(problem,
                               cons_val,
                               residuals,
                               working_set,
                               false);
}

SLEQP_RETCODE
sleqp_signed_feasibility_residuals(SleqpProblem* problem,
                                   const SleqpVec* cons_val,
                                   SleqpVec* residuals,
                                   SleqpWorkingSet* working_set)
{
  return feasibility_residuals(problem, cons_val, residuals, working_set, true);
}

SLEQP_RETCODE
sleqp_violation_inf_norm(SleqpProblem* problem,
                         SleqpVec* cons_val,
                         double* max_violation)
{
  const int num_constraints = sleqp_problem_num_cons(problem);

  const SleqpVec* lb = sleqp_problem_cons_lb(problem);
  const SleqpVec* ub = sleqp_problem_cons_ub(problem);
  SleqpVec* c        = cons_val;

  const int dim = c->dim;

  assert(dim == num_constraints);
  assert(dim == lb->dim);
  assert(dim == ub->dim);

  (*max_violation) = 0.;

  int k_c = 0, k_lb = 0, k_ub = 0;

  while (k_c < c->nnz || k_lb < lb->nnz || k_ub < ub->nnz)
  {
    double c_val = 0., lb_val = 0., ub_val = 0.;

    bool valid_c  = (k_c < c->nnz);
    bool valid_lb = (k_lb < lb->nnz);
    bool valid_ub = (k_ub < ub->nnz);

    int idx = valid_c ? c->indices[k_c] : dim + 1;
    idx     = SLEQP_MIN(idx, valid_lb ? lb->indices[k_lb] : dim + 1);
    idx     = SLEQP_MIN(idx, valid_ub ? ub->indices[k_ub] : dim + 1);

    if (valid_c && idx == c->indices[k_c])
    {
      c_val = c->data[k_c++];
    }

    if (valid_lb && idx == lb->indices[k_lb])
    {
      lb_val = lb->data[k_lb++];
    }

    if (valid_ub && idx == ub->indices[k_ub])
    {
      ub_val = ub->data[k_ub++];
    }

    const double upper_violation = SLEQP_MAX(c_val - ub_val, 0.);
    const double lower_violation = SLEQP_MAX(lb_val - c_val, 0.);

    if (upper_violation > 0.)
    {
      (*max_violation) = SLEQP_MAX(upper_violation, (*max_violation));
    }
    else if (lower_violation > 0.)
    {
      (*max_violation) = SLEQP_MAX(lower_violation, (*max_violation));
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_violation_one_norm(SleqpProblem* problem,
                         SleqpVec* cons_val,
                         double* total_violation)
{
  const int num_constraints = sleqp_problem_num_cons(problem);

  const SleqpVec* lb = sleqp_problem_cons_lb(problem);
  const SleqpVec* ub = sleqp_problem_cons_ub(problem);
  SleqpVec* c        = cons_val;

  const int dim = c->dim;

  assert(dim == num_constraints);
  assert(dim == lb->dim);
  assert(dim == ub->dim);

  (*total_violation) = 0.;

  int k_c = 0, k_lb = 0, k_ub = 0;

  while (k_c < c->nnz || k_lb < lb->nnz || k_ub < ub->nnz)
  {
    double c_val = 0., lb_val = 0., ub_val = 0.;

    bool valid_c  = (k_c < c->nnz);
    bool valid_lb = (k_lb < lb->nnz);
    bool valid_ub = (k_ub < ub->nnz);

    int idx = valid_c ? c->indices[k_c] : dim + 1;
    idx     = SLEQP_MIN(idx, valid_lb ? lb->indices[k_lb] : dim + 1);
    idx     = SLEQP_MIN(idx, valid_ub ? ub->indices[k_ub] : dim + 1);

    if (valid_c && idx == c->indices[k_c])
    {
      c_val = c->data[k_c++];
    }

    if (valid_lb && idx == lb->indices[k_lb])
    {
      lb_val = lb->data[k_lb++];
    }

    if (valid_ub && idx == ub->indices[k_ub])
    {
      ub_val = ub->data[k_ub++];
    }

    const double upper_violation = SLEQP_MAX(c_val - ub_val, 0.);
    const double lower_violation = SLEQP_MAX(lb_val - c_val, 0.);

    if (upper_violation > 0.)
    {
      (*total_violation) += upper_violation;
    }
    else if (lower_violation > 0.)
    {
      (*total_violation) += lower_violation;
    }
  }

  return SLEQP_OKAY;
}
