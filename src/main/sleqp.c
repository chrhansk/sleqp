#include "sleqp.h"

#include <assert.h>

#include "sleqp_cmp.h"
#include "sleqp_mem.h"

#include "lp/sleqp_lpi.h"
#include "lp/sleqp_lpi_soplex.h"

struct SleqpSolver
{
  SleqpProblem* problem;

  SleqpSparseVec* x;

  SleqpLPi* lp_interface;

  SleqpSparseVec* func_grad;
  SleqpSparseMatrix* cons_jac;
};

// copy and adjust upper / lower bounds
static SLEQP_RETCODE adjust_bounds(SleqpSparseVec* lb,
                                   SleqpSparseVec* ub,
                                   SleqpSparseVec** adj_lbstar,
                                   SleqpSparseVec** adj_ubstar)
{
  assert(lb->dim == ub->dim);

  const size_t dim = lb->dim;

  SLEQP_CALL(sleqp_sparse_vector_create(adj_lbstar,
                                        lb->dim,
                                        lb->nnz));

  SLEQP_CALL(sleqp_sparse_vector_create(adj_ubstar,
                                        ub->dim,
                                        ub->nnz));

  SleqpSparseVec* adj_lb = *adj_lbstar;
  SleqpSparseVec* adj_ub = *adj_ubstar;

  size_t k_lb = 0, k_ub = 0;

  while(k_lb < lb->nnz || k_ub < ub->nnz)
  {
    double lb_val = 0;
    double ub_val = 0;

    SLEQP_Bool valid_lb = (k_lb < lb->nnz);
    SLEQP_Bool valid_ub = (k_ub < ub->nnz);

    size_t lb_i = valid_lb ? lb->indices[k_lb] : lb->dim + 1;
    size_t ub_i = valid_ub ? ub->indices[k_ub] : ub->dim + 1;

    size_t idx = SLEQP_MIN(lb_i, ub_i);

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
  sleqp_malloc(star);

  SleqpProblem* sleqp = *star;

  assert(var_lb->dim == var_ub->dim);
  assert(cons_lb->dim == cons_ub->dim);

  sleqp->func = func;

  SleqpSparseVec** adj_lbstar;
  SleqpSparseVec** adj_ubstar;

  SLEQP_CALL(adjust_bounds(var_lb,
                           var_ub,
                           adj_lbstar,
                           adj_ubstar));

  sleqp->var_lb = *adj_lbstar;
  sleqp->var_ub = *adj_ubstar;

  SLEQP_CALL(adjust_bounds(cons_lb,
                           cons_ub,
                           adj_lbstar,
                           adj_ubstar));

  sleqp->cons_lb = *adj_lbstar;
  sleqp->cons_ub = *adj_ubstar;

  sleqp->num_variables = var_lb->dim;
  sleqp->num_constraints = cons_lb->dim;

  return SLEQP_OKAY;
}

#if 0
static SLEQP_RETCODE iterate(Sleqp* sleqp)
{
  size_t func_grad_nnz = 0;
  size_t cons_jac_nnz = 0;

  SLEQP_CALL(sleqp_func_set_value(sleqp->func,
                                  sleqp->x,
                                  sleqp->num_variables,
                                  &func_grad_nnz,
                                  &cons_jac_nnz));

  SLEQP_CALL(sleqp_sparse_vector_reserve(sleqp->func_grad,
                                         func_grad_nnz));

  /*
   * Reserve a litte more so we can add the two
   * identity matrices afterwards
   */
  SLEQP_CALL(sleqp_sparse_matrix_reserve(sleqp->cons_jac,
                                         cons_jac_nnz + 2*sleqp->num_constraints));

  /*
  SLEQP_CALL(sleqp_func_eval(sleqp->func,
                             sleqp->num_variables,
                             NULL,
                             fvals,
                             grad));
  */

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_problem_solve(Sleqp* sleqp,
                                  SleqpSparseVec* x)
{
  SleqpSparseVec* xclip;

  SLEQP_CALL(sleqp_sparse_vector_clip(x, sleqp->var_lb, sleqp->var_ub, &xclip));

  sleqp->x = xclip;

  iterate(sleqp);

  return SLEQP_OKAY;
}

#endif

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


SLEQP_RETCODE sleqp_solver_create(SleqpSolver** star,
                                  SleqpProblem* problem)
{
  sleqp_malloc(star);

  SleqpSolver* solver = *star;

  solver->problem = problem;

  SLEQP_CALL(sleqp_sparse_vector_create(&solver->func_grad,
                                        problem->num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_matrix_create(&solver->cons_jac,
                                        problem->num_constraints,
                                        problem->num_variables,
                                        0));

  // TODO: make this generic at a later point...
  SLEQP_CALL(sleqp_lpi_soplex_create_interface(&solver->lp_interface));

  size_t num_lp_variables = problem->num_variables + 2* problem->num_constraints;

  SLEQP_CALL(sleqp_lpi_create_problem(solver->lp_interface,
                                      num_lp_variables,
                                      problem->num_constraints));

  solver->x = NULL;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solve(SleqpSolver* solver,
                          SleqpSparseVec* x)
{
  SleqpSparseVec* xclip;

  SLEQP_CALL(sleqp_sparse_vector_clip(x,
                                      solver->problem->var_lb,
                                      solver->problem->var_ub,
                                      &xclip));

  solver->x = xclip;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solver_free(SleqpSolver** star)
{
  SleqpSolver* solver = *star;

  SLEQP_CALL(sleqp_lpi_free(&solver->lp_interface));

  SLEQP_CALL(sleqp_sparse_matrix_free(&solver->cons_jac));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->func_grad));

  if(solver->x)
  {
    SLEQP_CALL(sleqp_sparse_vector_free(&solver->x));
  }

  sleqp_free(star);

  return SLEQP_OKAY;
}
