#include "sleqp_util.h"

#include <assert.h>

#include "sleqp_cmp.h"

SLEQP_RETCODE sleqp_set_and_evaluate(SleqpProblem* problem,
                                     SleqpIterate* iterate,
                                     SLEQP_VALUE_REASON reason)
{
  int func_grad_nnz = 0;
  int cons_val_nnz = 0;
  int cons_jac_nnz = 0;

  SLEQP_CALL(sleqp_func_set_value(problem->func,
                                  sleqp_iterate_get_primal(iterate),
                                  reason,
                                  &func_grad_nnz,
                                  &cons_val_nnz,
                                  &cons_jac_nnz));

  SleqpSparseVec* func_grad = sleqp_iterate_get_func_grad(iterate);
  SleqpSparseMatrix* cons_jac = sleqp_iterate_get_cons_jac(iterate);
  SleqpSparseVec* cons_val = sleqp_iterate_get_cons_val(iterate);

  SLEQP_CALL(sleqp_sparse_vector_reserve(func_grad,
                                         func_grad_nnz));

  SLEQP_CALL(sleqp_sparse_vector_reserve(cons_val,
                                         cons_val_nnz));

  SLEQP_CALL(sleqp_sparse_matrix_reserve(cons_jac,
                                         cons_jac_nnz));

  double func_val;

  SLEQP_CALL(sleqp_func_eval(problem->func,
                             NULL,
                             &func_val,
                             func_grad,
                             cons_val,
                             cons_jac));

  SLEQP_CALL(sleqp_iterate_set_func_val(iterate, func_val));

  assert(sleqp_sparse_vector_valid(func_grad));
  assert(sleqp_sparse_vector_valid(cons_val));

  assert(sleqp_sparse_matrix_valid(cons_jac));

  return SLEQP_OKAY;
}


SLEQP_RETCODE sleqp_get_violated_multipliers(SleqpProblem* problem,
                                             SleqpSparseVec* x,
                                             SleqpSparseVec* cons_vals,
                                             double penalty_parameter,
                                             SleqpSparseVec* multipliers,
                                             SleqpWorkingSet* working_set,
                                             double eps)
{
  SleqpSparseVec* lb = problem->cons_lb;
  SleqpSparseVec* ub = problem->cons_ub;
  SleqpSparseVec* v = cons_vals;

  SLEQP_CALL(sleqp_sparse_vector_clear(multipliers));

  // TODO: use active set cons size instead...
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

    if(sleqp_gt(c_val, ub_val, eps))
    {
      SLEQP_CALL(sleqp_sparse_vector_push(multipliers,
                                          idx,
                                          penalty_parameter));
    }
    else if(sleqp_lt(c_val, lb_val, eps))
    {
      SLEQP_CALL(sleqp_sparse_vector_push(multipliers,
                                          idx,
                                          -penalty_parameter));
    }

  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_max_step_length(SleqpSparseVec* x,
                                    SleqpSparseVec* d,
                                    SleqpSparseVec* l,
                                    SleqpSparseVec* u,
                                    double* max_step_length)
{
  const int dim = x->dim;

  assert(l->dim == dim);
  assert(u->dim == dim);
  assert(d->dim == dim);

  assert((*max_step_length) > 0.);

  // upper bound

  {
    int k_x = 0, k_d = 0, k_u = 0;

    while(k_x < x->nnz || k_d < d->nnz || k_u < u->nnz)
    {
      bool valid_x = (k_x < x->nnz);
      bool valid_d = (k_d < d->nnz);
      bool valid_u = (k_u < u->nnz);

      int i_combined = valid_x ? x->indices[k_x] : dim + 1;
      i_combined = SLEQP_MIN(i_combined, valid_d ? d->indices[k_d] : dim + 1);
      i_combined = SLEQP_MIN(i_combined, valid_u ? u->indices[k_u] : dim + 1);

      valid_x = valid_x && x->indices[k_x] == i_combined;
      valid_d = valid_d && d->indices[k_d] == i_combined;
      valid_u = valid_u && u->indices[k_u] == i_combined;

      double x_value = valid_x ? x->data[k_x] : 0.;
      double d_value = valid_d ? d->data[k_d] : 0.;
      double u_value = valid_u ? u->data[k_u] : 0.;

      double diff = u_value - x_value;

      if((d_value > 0. && diff > 0.) ||
         (d_value < 0. && diff < 0.))
      {
        double current_bound = diff / d_value;

        (*max_step_length) = SLEQP_MIN(*max_step_length, current_bound);
      }

      if(valid_x)
      {
        ++k_x;
      }

      if(valid_d)
      {
        ++k_d;
      }

      if(valid_u)
      {
        ++k_u;
      }
    }
  }

  // lower bound

  {
    int k_x = 0, k_d = 0, k_l = 0;

    while(k_x < x->nnz || k_d < d->nnz || k_l < l->nnz)
    {
      bool valid_x = (k_x < x->nnz);
      bool valid_d = (k_d < d->nnz);
      bool valid_l = (k_l < l->nnz);

      int i_combined = valid_x ? x->indices[k_x] : dim + 1;
      i_combined = SLEQP_MIN(i_combined, valid_d ? d->indices[k_d] : dim + 1);
      i_combined = SLEQP_MIN(i_combined, valid_l ? l->indices[k_l] : dim + 1);

      valid_x = valid_x && x->indices[k_x] == i_combined;
      valid_d = valid_d && d->indices[k_d] == i_combined;
      valid_l = valid_l && l->indices[k_l] == i_combined;

      double x_value = valid_x ? x->data[k_x] : 0.;
      double d_value = valid_d ? d->data[k_d] : 0.;
      double l_value = valid_l ? l->data[k_l] : 0.;

      double diff = l_value - x_value;

      if((d_value < 0. && diff < 0.) ||
         (d_value > 0. && diff > 0.))
      {
        double current_bound = diff / d_value;

        (*max_step_length) = SLEQP_MIN(*max_step_length, current_bound);
      }

      if(valid_x)
      {
        ++k_x;
      }

      if(valid_d)
      {
        ++k_d;
      }

      if(valid_l)
      {
        ++k_l;
      }
    }
  }

  assert((*max_step_length) >= 0.);

  return SLEQP_OKAY;
}


SLEQP_RETCODE sleqp_get_violation(SleqpProblem* problem,
                                  SleqpIterate* iterate,
                                  double eps,
                                  SleqpSparseVec* violation)
{
  int num_constraints = problem->num_constraints;

  SleqpSparseVec* lb = problem->cons_lb;
  SleqpSparseVec* ub = problem->cons_ub;
  SleqpSparseVec* v = sleqp_iterate_get_cons_val(iterate);

  const int dim = v->dim;

  assert(dim == num_constraints);
  assert(dim == lb->dim);
  assert(dim == ub->dim);
  assert(dim == violation->dim);

  {
    int max_size = lb->nnz + ub->nnz + v->nnz;
    max_size = SLEQP_MIN(max_size, num_constraints);

    SLEQP_CALL(sleqp_sparse_vector_clear(violation));
    SLEQP_CALL(sleqp_sparse_vector_reserve(violation, max_size));
  }

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

    if(sleqp_gt(c_val, ub_val, eps))
    {
      SLEQP_CALL(sleqp_sparse_vector_push(violation,
                                          idx,
                                          c_val - ub_val));
    }
    else if(sleqp_lt(c_val, lb_val, eps))
    {
      SLEQP_CALL(sleqp_sparse_vector_push(violation,
                                          idx,
                                          c_val - lb_val));
    }

  }

  return SLEQP_OKAY;
}
