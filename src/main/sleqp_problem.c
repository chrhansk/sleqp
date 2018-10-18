#include "sleqp_problem.h"

#include <assert.h>

#include "sleqp_cmp.h"
#include "sleqp_mem.h"

static SLEQP_RETCODE adjust_bounds(SleqpSparseVec* lb,
                                   SleqpSparseVec* ub,
                                   SleqpSparseVec** adj_lbstar,
                                   SleqpSparseVec** adj_ubstar)
{
  assert(lb->dim == ub->dim);

  SLEQP_CALL(sleqp_sparse_vector_create(adj_lbstar,
                                        lb->dim,
                                        lb->nnz));

  SLEQP_CALL(sleqp_sparse_vector_create(adj_ubstar,
                                        ub->dim,
                                        ub->nnz));

  SleqpSparseVec* adj_lb = *adj_lbstar;
  SleqpSparseVec* adj_ub = *adj_ubstar;

  int k_lb = 0, k_ub = 0;

  while(k_lb < lb->nnz || k_ub < ub->nnz)
  {
    double lb_val = 0;
    double ub_val = 0;

    SLEQP_Bool valid_lb = (k_lb < lb->nnz);
    SLEQP_Bool valid_ub = (k_ub < ub->nnz);

    int lb_i = valid_lb ? lb->indices[k_lb] : lb->dim + 1;
    int ub_i = valid_ub ? ub->indices[k_ub] : ub->dim + 1;

    int idx = SLEQP_MIN(lb_i, ub_i);

    if(valid_lb && idx == lb_i)
    {
      lb_val = lb->data[k_lb];
    }

    if(valid_ub && idx == ub_i)
    {
      ub_val = ub->data[k_ub];
    }

    if(sleqp_gt(lb_val, ub_val))
    {
      SLEQP_CALL(sleqp_sparse_vector_free(adj_lbstar));

      SLEQP_CALL(sleqp_sparse_vector_free(adj_ubstar));

      return SLEQP_INVALID;
    }

    if(valid_lb)
    {
      SLEQP_CALL(sleqp_sparse_vector_push(adj_lb, k_lb, SLEQP_MIN(lb_val, ub_val)));

      ++k_lb;
    }

    if(valid_ub)
    {
      SLEQP_CALL(sleqp_sparse_vector_push(adj_ub, k_ub, SLEQP_MAX(lb_val, ub_val)));

      ++k_ub;
    }
  }


  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_problem_create(SleqpProblem** star,
                                   SleqpFunc* func,
                                   SleqpSparseVec* var_lb,
                                   SleqpSparseVec* var_ub,
                                   SleqpSparseVec* cons_lb,
                                   SleqpSparseVec* cons_ub)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpProblem* sleqp = *star;

  assert(var_lb->dim == var_ub->dim);
  assert(cons_lb->dim == cons_ub->dim);

  assert(sleqp_sparse_vector_valid(var_lb));
  assert(sleqp_sparse_vector_valid(var_ub));

  assert(sleqp_sparse_vector_valid(cons_lb));
  assert(sleqp_sparse_vector_valid(cons_ub));

  sleqp->func = func;

  SleqpSparseVec* adj_lb;
  SleqpSparseVec* adj_ub;

  SLEQP_CALL(adjust_bounds(var_lb,
                           var_ub,
                           &adj_lb,
                           &adj_ub));

  sleqp->var_lb = adj_lb;
  sleqp->var_ub = adj_ub;

  SLEQP_CALL(adjust_bounds(cons_lb,
                           cons_ub,
                           &adj_lb,
                           &adj_ub));

  sleqp->cons_lb = adj_lb;
  sleqp->cons_ub = adj_ub;

  sleqp->num_variables = var_lb->dim;
  sleqp->num_constraints = cons_lb->dim;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_problem_free(SleqpProblem** star)
{
  SleqpProblem* sleqp = *star;

  sleqp_sparse_vector_free(&sleqp->var_lb);
  sleqp_sparse_vector_free(&sleqp->var_ub);

  sleqp_sparse_vector_free(&sleqp->cons_lb);
  sleqp_sparse_vector_free(&sleqp->cons_ub);

  sleqp_free(star);

  return SLEQP_OKAY;
}
