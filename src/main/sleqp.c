#include "sleqp.h"

#include <assert.h>

#include "sleqp_cmp.h"
#include "sleqp_mem.h"

#include "sleqp_cauchy.h"
#include "sleqp_iterate.h"

#include "lp/sleqp_lpi.h"
#include "lp/sleqp_lpi_soplex.h"

struct SleqpSolver
{
  SleqpProblem* problem;

  SleqpIterate* iterate;

  SleqpLPi* lp_interface;

  SleqpCauchyData* cauchy_data;

  double trust_radius;
  double penalty;
};

// copy and adjust upper / lower bounds
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
  SLEQP_CALL(sleqp_malloc(star));

  SleqpProblem* sleqp = *star;

  assert(var_lb->dim == var_ub->dim);
  assert(cons_lb->dim == cons_ub->dim);

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

SLEQP_RETCODE sleqp_set_and_evaluate(SleqpProblem* problem,
                                     SleqpIterate* iterate)
{
  size_t func_grad_nnz = 0;
  size_t cons_val_nnz = 0;
  size_t cons_jac_nnz = 0;

  SLEQP_CALL(sleqp_func_set_value(problem->func,
                                  iterate->x,
                                  problem->num_variables,
                                  &func_grad_nnz,
                                  &cons_val_nnz,
                                  &cons_jac_nnz));

  /*
   * Reserve a litte more so we can add slack penalties
   * without reallocations
   */
  SLEQP_CALL(sleqp_sparse_vector_reserve(iterate->func_grad,
                                         func_grad_nnz + 2*problem->num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_reserve(iterate->cons_val,
                                         cons_val_nnz));

  /*
   * Reserve a litte more so we can add the two
   * identity matrices afterwards
   */
  SLEQP_CALL(sleqp_sparse_matrix_reserve(iterate->cons_jac,
                                         cons_jac_nnz + 2*problem->num_constraints));

  SLEQP_CALL(sleqp_func_eval(problem->func,
                             problem->num_variables,
                             NULL,
                             &iterate->func_val,
                             iterate->func_grad,
                             iterate->cons_val,
                             iterate->cons_jac));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE iterate(SleqpSolver* solver)
{
  /*

  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;



  SLEQP_CALL(sleqp_cauchy_compute_direction(problem,
                                            iterate,
                                            solver->cauchy_data,
                                            solver->penalty,
                                            solver->trust_radius));
  */

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


SLEQP_RETCODE sleqp_solver_create(SleqpSolver** star,
                                  SleqpProblem* problem,
                                  SleqpSparseVec* x)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpSolver* solver = *star;
  size_t num_constraints = problem->num_constraints;
  size_t num_variables = problem->num_variables;

  solver->problem = problem;

  SleqpSparseVec* xclip;

  SLEQP_CALL(sleqp_sparse_vector_clip(x,
                                      solver->problem->var_lb,
                                      solver->problem->var_ub,
                                      &xclip));

  SLEQP_CALL(sleqp_iterate_create(&solver->iterate,
                                  solver->problem,
                                  xclip));

  // TODO: make this generic at a later point...
  SLEQP_CALL(sleqp_lpi_soplex_create_interface(&solver->lp_interface));

  SLEQP_CALL(sleqp_lpi_create_problem(solver->lp_interface,
                                      num_variables + 2*num_constraints,
                                      num_constraints));

  SLEQP_CALL(sleqp_cauchy_data_create(&solver->cauchy_data,
                                      problem,
                                      solver->lp_interface));

  // TODO: Set this to something different?!
  solver->trust_radius = 1.;
  solver->penalty = 1.;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solve(SleqpSolver* solver)
{
  SLEQP_CALL(iterate(solver));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solver_free(SleqpSolver** star)
{
  SleqpSolver* solver = *star;

  SLEQP_CALL(sleqp_cauchy_data_free(&solver->cauchy_data));

  SLEQP_CALL(sleqp_lpi_free(&solver->lp_interface));

  SLEQP_CALL(sleqp_iterate_free(&solver->iterate));

  sleqp_free(star);

  return SLEQP_OKAY;
}
