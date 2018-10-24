#include "sleqp.h"

#include <assert.h>

#include "sleqp_cmp.h"

SLEQP_RETCODE sleqp_set_and_evaluate(SleqpProblem* problem,
                                     SleqpIterate* iterate)
{
  int func_grad_nnz = 0;
  int cons_val_nnz = 0;
  int cons_jac_nnz = 0;

  SLEQP_CALL(sleqp_func_set_value(problem->func,
                                  iterate->x,
                                  &func_grad_nnz,
                                  &cons_val_nnz,
                                  &cons_jac_nnz));

  SLEQP_CALL(sleqp_sparse_vector_reserve(iterate->func_grad, func_grad_nnz));

  SLEQP_CALL(sleqp_sparse_vector_reserve(iterate->cons_val, cons_val_nnz));

  SLEQP_CALL(sleqp_sparse_matrix_reserve(iterate->cons_jac, cons_jac_nnz));

  SLEQP_CALL(sleqp_func_eval(problem->func,
                             NULL,
                             &iterate->func_val,
                             iterate->func_grad,
                             iterate->cons_val,
                             iterate->cons_jac));

  assert(sleqp_sparse_vector_valid(iterate->func_grad));
  assert(sleqp_sparse_vector_valid(iterate->cons_val));

  assert(sleqp_sparse_matrix_valid(iterate->cons_jac));

  return SLEQP_OKAY;
}


SLEQP_RETCODE sleqp_get_violated_constraints(SleqpProblem* problem,
                                             SleqpSparseVec* x,
                                             SleqpSparseVec* cons_vals,
                                             double penalty_parameter,
                                             SleqpSparseVec* multipliers,
                                             SleqpActiveSet* active_set)
{
  SleqpSparseVec* lb = problem->cons_lb;
  SleqpSparseVec* ub = problem->cons_ub;
  SleqpSparseVec* v = cons_vals;

  SLEQP_ACTIVE_STATE* cons_states = active_set ? sleqp_active_set_cons_states(active_set) : NULL;

  const int dim = v->dim;

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

    if(cons_states && cons_states[idx] != SLEQP_INACTIVE)
    {
      continue;
    }

    if(sleqp_gt(c_val, ub_val))
    {
      SLEQP_CALL(sleqp_sparse_vector_push(multipliers,
                                          idx,
                                          penalty_parameter));
    }
    else if(sleqp_lt(c_val, lb_val))
    {
      SLEQP_CALL(sleqp_sparse_vector_push(multipliers,
                                          idx,
                                          -penalty_parameter));
    }

  }

  return SLEQP_OKAY;
}
