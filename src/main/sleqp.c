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


SLEQP_RETCODE sleqp_set_and_evaluate(SleqpProblem* problem,
                                     SleqpIterate* iterate)
{
  size_t func_grad_nnz = 0;
  size_t cons_val_nnz = 0;
  size_t cons_jac_nnz = 0;

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
