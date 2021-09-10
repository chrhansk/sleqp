#include "lsqr.h"

#include <math.h>

#include "aug_jacobian.h"
#include "cmp.h"
#include "fail.h"
#include "feas.h"
#include "log.h"
#include "lsq.h"
#include "mem.h"
#include "util.h"
#include "working_step.h"

#include "tr/tr_util.h"

static const double tolerance_factor = 1e-2;

typedef struct {
  SleqpSparseVec* u;
  SleqpSparseVec* v;
  SleqpSparseVec* w;

  SleqpSparseVec* p;
  SleqpSparseVec* q;

  SleqpSparseVec* x;

  SleqpSparseVec* d;
} LSQR;

typedef struct {
  SleqpSparseVec* projected_direction;
  SleqpSparseVec* scaled_cons_product;
  SleqpSparseVec* lsq_product;
} Forward;

typedef struct {
  SleqpSparseVec* product;

  SleqpSparseVec* lsq_direction;
  SleqpSparseVec* cons_direction;

  SleqpSparseVec* lsq_product;
  SleqpSparseVec* cons_product;

} Adjoint;

struct SleqpLSQRSolver {
  int refcount;

  SleqpProblem* problem;
  SleqpWorkingStep* working_step;
  SleqpParams* params;

  SleqpIterate* iterate;

  SleqpSparseVec* lsq_residuals;
  // SleqpSparseVec* lsq_step_forward;
  SleqpSparseVec* lsq_rhs;
  SleqpSparseVec* rhs;

  SleqpSparseVec* linear_cons_val;
  SleqpSparseVec* linear_cons_residuals;

  SleqpSparseMatrix* scaled_violated_cons_jac;
  SleqpSparseVec* scaled_cons_residuals;
  // SleqpSparseVec* scaled_cons_product;

  int* removed_cons;
  int num_removed_cons;

  int forward_dim;
  int adjoint_dim;
  int num_violated_cons;

  SleqpSparseVec* sparse_cache;
  double* dense_cache;

  // SleqpSparseVec* projected_direction;

  SleqpAugJacobian* jacobian;

  double trust_radius;
  double penalty_parameter;

  Forward forward;
  Adjoint adjoint;
  LSQR lsqr;
};

SLEQP_RETCODE sleqp_lsqr_solver_create(SleqpLSQRSolver** star,
                                       SleqpProblem* problem,
                                       SleqpWorkingStep* step,
                                       SleqpParams* params)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpLSQRSolver* solver = *star;

  *solver = (SleqpLSQRSolver) {0};

  solver->refcount = 1;

  SLEQP_CALL(sleqp_problem_capture(problem));
  solver->problem = problem;

  SleqpFunc* func = sleqp_problem_func(problem);

  assert(sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_LSQ);

  SLEQP_CALL(sleqp_working_step_capture(step));
  solver->working_step = step;

  SLEQP_CALL(sleqp_params_capture(params));
  solver->params = params;

  const int num_variables = sleqp_problem_num_variables(problem);
  const int num_constraints = sleqp_problem_num_constraints(problem);
  const int num_residuals = sleqp_lsq_func_num_residuals(func);

  solver->forward_dim = num_variables;
  solver->adjoint_dim = 0;

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->lsq_residuals,
                                              num_residuals));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->lsq_rhs,
                                              num_residuals));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->rhs,
                                              0));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->linear_cons_val,
                                              num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->linear_cons_residuals,
                                              num_constraints));

  SLEQP_CALL(sleqp_sparse_matrix_create(&solver->scaled_violated_cons_jac,
                                        0,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->scaled_cons_residuals,
                                              0));

  SLEQP_CALL(sleqp_alloc_array(&solver->removed_cons,
                               num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->sparse_cache,
                                              num_constraints));

  SLEQP_CALL(sleqp_alloc_array(&solver->dense_cache,
                               num_constraints));

  {
    SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->lsqr.u,
                                                0));

    SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->lsqr.v,
                                                solver->forward_dim));

    SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->lsqr.w,
                                                solver->forward_dim));

    SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->lsqr.p,
                                                0));

    SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->lsqr.q,
                                                solver->forward_dim));

    SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->lsqr.x,
                                                solver->forward_dim));

    SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->lsqr.d,
                                                solver->forward_dim));
  }

  {
    SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->forward.lsq_product,
                                                num_residuals));

    SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->forward.scaled_cons_product,
                                                0));

    SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->forward.projected_direction,
                                                num_variables));
  }

  {
    SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->adjoint.lsq_direction,
                                                num_residuals));

    SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->adjoint.cons_direction,
                                                0));

    SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->adjoint.lsq_product,
                                                solver->forward_dim));

    SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->adjoint.cons_product,
                                                solver->forward_dim));

    SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->adjoint.product,
                                                solver->forward_dim));
  }

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE compute_lsq_rhs(SleqpLSQRSolver* solver)
{
  SleqpProblem* problem = solver->problem;
  SleqpFunc* func = sleqp_problem_func(problem);

  const double zero_eps = sleqp_params_get(solver->params,
                                           SLEQP_PARAM_ZERO_EPS);

  SleqpSparseVec* initial_step = sleqp_working_step_get_step(solver->working_step);

  SLEQP_CALL(sleqp_lsq_func_residuals(func,
                                      solver->lsq_residuals));

  SLEQP_CALL(sleqp_lsq_func_jac_forward(func,
                                        initial_step,
                                        solver->forward.lsq_product));

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(solver->lsq_residuals,
                                            solver->forward.lsq_product,
                                            -1.,
                                            -1.,
                                            zero_eps,
                                            solver->lsq_rhs));

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE compute_cons_rhs(SleqpLSQRSolver* solver)
{
  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;

  const double zero_eps = sleqp_params_get(solver->params,
                                           SLEQP_PARAM_ZERO_EPS);

  const int num_constraints = sleqp_problem_num_constraints(problem);

  SleqpSparseVec* cons_val = sleqp_iterate_get_cons_val(iterate);
  SleqpSparseMatrix* cons_jac = sleqp_iterate_get_cons_jac(iterate);
  SleqpWorkingSet* working_set = sleqp_iterate_get_working_set(iterate);

  SleqpSparseVec* initial_step = sleqp_working_step_get_step(solver->working_step);

  SLEQP_CALL(sleqp_sparse_matrix_vector_product(cons_jac,
                                                initial_step,
                                                solver->dense_cache));

  SLEQP_CALL(sleqp_sparse_vector_from_raw(solver->sparse_cache,
                                          solver->dense_cache,
                                          num_constraints,
                                          zero_eps));

  SLEQP_CALL(sleqp_sparse_vector_add(cons_val,
                                     solver->sparse_cache,
                                     zero_eps,
                                     solver->linear_cons_val));

  SLEQP_CALL(sleqp_feasibility_residuals(problem,
                                         solver->linear_cons_val,
                                         solver->linear_cons_residuals,
                                         working_set));

  SLEQP_CALL(sleqp_sparse_vector_clear(solver->scaled_cons_residuals));

  SLEQP_CALL(sleqp_sparse_vector_resize(solver->scaled_cons_residuals,
                                        solver->num_violated_cons));

  SLEQP_CALL(sleqp_sparse_vector_remove_entries(solver->linear_cons_residuals,
                                                solver->scaled_cons_residuals,
                                                solver->removed_cons,
                                                solver->num_removed_cons));

  SLEQP_CALL(sleqp_sparse_vector_scale(solver->scaled_cons_residuals,
                                       solver->penalty_parameter));

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE compute_cons_matrix(SleqpLSQRSolver* solver)
{
  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;

  const int num_variables = sleqp_problem_num_variables(problem);

  SleqpSparseMatrix* cons_jac = sleqp_iterate_get_cons_jac(iterate);

  SLEQP_CALL(sleqp_sparse_matrix_clear(solver->scaled_violated_cons_jac));

  SLEQP_CALL(sleqp_sparse_matrix_resize(solver->scaled_violated_cons_jac,
                                        solver->num_violated_cons,
                                        num_variables));

  SLEQP_CALL(sleqp_sparse_matrix_remove_rows(cons_jac,
                                             solver->scaled_violated_cons_jac,
                                             solver->removed_cons,
                                             solver->num_removed_cons));

  SLEQP_CALL(sleqp_sparse_matrix_scale(solver->scaled_violated_cons_jac,
                                       solver->penalty_parameter));

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE compute_removed_cons(SleqpLSQRSolver* solver)
{
  SleqpProblem* problem = solver->problem;
  const int num_constraints = sleqp_problem_num_constraints(problem);

  SleqpSparseVec* violated_cons_mult = sleqp_working_step_get_violated_cons_multipliers(solver->working_step);

  solver->num_removed_cons = 0;

  int k = 0;

  for(int i = 0; i < num_constraints; ++i)
  {
    if(k < violated_cons_mult->nnz && i == violated_cons_mult->indices[k])
    {
      ++k;
    }
    else
    {
      solver->removed_cons[solver->num_removed_cons++] = i;
    }
  }

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE resize_to_dimensions(SleqpLSQRSolver* solver)
{
  SLEQP_CALL(sleqp_sparse_vector_resize(solver->adjoint.cons_direction,
                                        solver->num_violated_cons));

  SLEQP_CALL(sleqp_sparse_vector_resize(solver->lsqr.u,
                                        solver->adjoint_dim));

  SLEQP_CALL(sleqp_sparse_vector_resize(solver->lsqr.p,
                                        solver->adjoint_dim));

  SLEQP_CALL(sleqp_sparse_vector_resize(solver->rhs,
                                        solver->adjoint_dim));

  SLEQP_CALL(sleqp_sparse_vector_resize(solver->forward.scaled_cons_product,
                                        solver->num_violated_cons));

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE set_adjoint_dimension(SleqpLSQRSolver* solver)
{
  SleqpProblem* problem = solver->problem;
  SleqpFunc* func = sleqp_problem_func(problem);

  const int num_constraints = sleqp_problem_num_constraints(problem);
  const int num_residuals = sleqp_lsq_func_num_residuals(func);

  solver->num_violated_cons = num_constraints - solver->num_removed_cons;
  solver->adjoint_dim = num_residuals + solver->num_violated_cons;

  SLEQP_CALL(resize_to_dimensions(solver));

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE compute_rhs(SleqpLSQRSolver* solver)
{
  SLEQP_CALL(compute_lsq_rhs(solver));

  SLEQP_CALL(compute_cons_rhs(solver));

  SLEQP_CALL(sleqp_sparse_vector_clear(solver->rhs));

  SLEQP_CALL(sleqp_sparse_vector_concat(solver->lsq_rhs,
                                        solver->scaled_cons_residuals,
                                        solver->rhs));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_lsqr_solver_set_iterate(SleqpLSQRSolver* solver,
                                            SleqpIterate* iterate,
                                            SleqpAugJacobian* jacobian,
                                            double trust_radius,
                                            double penalty_parameter)
{
  SLEQP_CALL(sleqp_working_step_set_iterate(solver->working_step,
                                            iterate,
                                            jacobian,
                                            trust_radius));

  {
    SLEQP_CALL(sleqp_iterate_release(&solver->iterate));

    SLEQP_CALL(sleqp_iterate_capture(iterate));

    solver->iterate = iterate;
  }

  {
    SLEQP_CALL(sleqp_aug_jacobian_release(&solver->jacobian));

    SLEQP_CALL(sleqp_aug_jacobian_capture(jacobian));

    solver->jacobian = jacobian;
  }

  solver->trust_radius = sleqp_working_step_get_reduced_trust_radius(solver->working_step);
  solver->penalty_parameter = penalty_parameter;

  SLEQP_CALL(compute_removed_cons(solver));

  SLEQP_CALL(set_adjoint_dimension(solver));

  SLEQP_CALL(compute_rhs(solver));

  SLEQP_CALL(compute_cons_matrix(solver));

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE forward_product(SleqpLSQRSolver* solver,
                              SleqpSparseVec* direction,
                              SleqpSparseVec* product)
{
  SleqpProblem* problem = solver->problem;
  SleqpFunc* func = sleqp_problem_func(problem);

  assert(direction->dim == solver->forward_dim);
  assert(product->dim == solver->adjoint_dim);

  const double zero_eps = sleqp_params_get(solver->params,
                                           SLEQP_PARAM_ZERO_EPS);

  SLEQP_CALL(sleqp_aug_jacobian_projection(solver->jacobian,
                                           direction,
                                           solver->forward.projected_direction,
                                           NULL));

  SLEQP_CALL(sleqp_lsq_func_jac_forward(func,
                                        solver->forward.projected_direction,
                                        solver->forward.lsq_product));

  SLEQP_CALL(sleqp_sparse_matrix_vector_product(solver->scaled_violated_cons_jac,
                                                solver->forward.projected_direction,
                                                solver->dense_cache));

  SLEQP_CALL(sleqp_sparse_vector_from_raw(solver->forward.scaled_cons_product,
                                          solver->dense_cache,
                                          solver->num_violated_cons,
                                          zero_eps));

  SLEQP_CALL(sleqp_sparse_vector_concat(solver->forward.lsq_product,
                                        solver->forward.scaled_cons_product,
                                        product));

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE adjoint_product(SleqpLSQRSolver* solver,
                              SleqpSparseVec* direction,
                              SleqpSparseVec* product)
{
  SleqpProblem* problem = solver->problem;
  SleqpFunc* func = sleqp_problem_func(problem);

  const double zero_eps = sleqp_params_get(solver->params,
                                           SLEQP_PARAM_ZERO_EPS);

  const int num_residuals = sleqp_lsq_func_num_residuals(func);

  assert(product->dim == solver->forward_dim);
  assert(direction->dim == solver->adjoint_dim);

  // split product
  {
    SLEQP_CALL(sleqp_sparse_vector_clear(solver->adjoint.lsq_direction));
    SLEQP_CALL(sleqp_sparse_vector_reserve(solver->adjoint.lsq_direction,
                                           SLEQP_MIN(num_residuals, direction->nnz)));

    SLEQP_CALL(sleqp_sparse_vector_clear(solver->adjoint.cons_direction));
    SLEQP_CALL(sleqp_sparse_vector_reserve(solver->adjoint.cons_direction,
                                           SLEQP_MIN(solver->num_violated_cons, direction->nnz)));

    for(int k = 0; k < direction->nnz; ++k)
    {
      const int i = direction->indices[k];
      const double v = direction->data[k];

      if(i < num_residuals)
      {
        SLEQP_CALL(sleqp_sparse_vector_push(solver->adjoint.lsq_direction,
                                            i,
                                            v));
      }
      else
      {
        SLEQP_CALL(sleqp_sparse_vector_push(solver->adjoint.cons_direction,
                                            i - num_residuals,
                                            v));
      }
    }
  }

  SLEQP_CALL(sleqp_lsq_func_jac_adjoint(func,
                                        solver->adjoint.lsq_direction,
                                        solver->adjoint.lsq_product));

  SLEQP_CALL(sleqp_sparse_matrix_trans_vector_product(solver->scaled_violated_cons_jac,
                                                      solver->adjoint.cons_direction,
                                                      zero_eps,
                                                      solver->adjoint.cons_product));

  SLEQP_CALL(sleqp_sparse_vector_add(solver->adjoint.lsq_product,
                                     solver->adjoint.cons_product,
                                     zero_eps,
                                     solver->adjoint.product));

  SLEQP_CALL(sleqp_aug_jacobian_projection(solver->jacobian,
                                           solver->adjoint.product,
                                           product,
                                           NULL));

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE normalize(SleqpSparseVec* vec,
                        double* norm)
{
  *norm = sleqp_sparse_vector_norm(vec);

  if(*norm)
  {
    SLEQP_CALL(sleqp_sparse_vector_scale(vec, 1./(*norm)));
  }

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE lsqr_loop(SleqpLSQRSolver* solver)
{
  SleqpProblem* problem = solver->problem;
  const int num_variables = sleqp_problem_num_variables(problem);
  int iteration = 0;

  const double eps = sleqp_params_get(solver->params,
                                      SLEQP_PARAM_EPS);

  const double zero_eps = sleqp_params_get(solver->params,
                                           SLEQP_PARAM_ZERO_EPS);

  const double stat_eps = sleqp_params_get(solver->params,
                                           SLEQP_PARAM_STATIONARITY_TOL);

  const double rel_tol = stat_eps * tolerance_factor;

  {
    sleqp_log_debug("Solving a LSQR subproblem with %d constraint residuals",
                    solver->num_violated_cons);
  }

  SleqpSparseVec* x = solver->lsqr.x;
  SleqpSparseVec* u = solver->lsqr.u;
  SleqpSparseVec* v = solver->lsqr.v;
  SleqpSparseVec* w = solver->lsqr.w;

  SleqpSparseVec* p = solver->lsqr.p;
  SleqpSparseVec* q = solver->lsqr.q;

  SleqpSparseVec* d = solver->lsqr.d;

  SleqpSparseVec* b = solver->rhs;
  SleqpSparseVec* t = solver->sparse_cache;

  double alpha, beta;

  SLEQP_CALL(sleqp_sparse_vector_copy(b, u));
  SLEQP_CALL(normalize(u, &beta));

  SLEQP_CALL(adjoint_product(solver, u, v));
  SLEQP_CALL(normalize(v, &alpha));

  SLEQP_CALL(sleqp_sparse_vector_copy(v, w));

  SLEQP_CALL(sleqp_sparse_vector_clear(x));

  double phib = beta;
  double rhob = alpha;

  for(iteration = 1; iteration <= num_variables; ++iteration)
  {
    // Continue bidiagonalization
    SLEQP_CALL(forward_product(solver, v, p));
    SLEQP_CALL(sleqp_sparse_vector_add_scaled(p, u, 1., -alpha, zero_eps, t));
    SLEQP_CALL(sleqp_sparse_vector_copy(t, u));
    SLEQP_CALL(normalize(u, &beta));

    SLEQP_CALL(adjoint_product(solver, u, q));
    SLEQP_CALL(sleqp_sparse_vector_add_scaled(q, v, 1., -beta, zero_eps, t));
    SLEQP_CALL(sleqp_sparse_vector_copy(t, v));
    SLEQP_CALL(normalize(v, &alpha));

    // Construct next orthogonal transformation (Givens rotation)
    double rho = hypot(rhob, beta);
    const double c = rhob / rho;
    const double s = beta / rho;
    double theta = s * alpha;
    rhob = (-c) * alpha;
    double phi = c * phib;
    phib = s * phib;

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(x,
                                              w,
                                              1.,
                                              phi / rho,
                                              zero_eps,
                                              t));

    const double norm = sleqp_sparse_vector_norm(t);

    if(sleqp_is_gt(norm, solver->trust_radius, eps))
    {
      SLEQP_CALL(sleqp_sparse_vector_add_scaled(x,
                                                t,
                                                -1.,
                                                1.,
                                                zero_eps,
                                                d));

      SLEQP_CALL(sleqp_tr_compute_bdry_sol(x,
                                           d,
                                           solver->params,
                                           solver->trust_radius,
                                           t));

      SLEQP_CALL(sleqp_sparse_vector_copy(t, x));

      sleqp_log_debug("LSQR solver terminated with a boundary solution after %d steps",
                      iteration);

      return SLEQP_OKAY;
    }

    SLEQP_CALL(sleqp_sparse_vector_copy(t, x));

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(v,
                                              w,
                                              1.,
                                              - theta / rho,
                                              zero_eps,
                                              t));

    SLEQP_CALL(sleqp_sparse_vector_copy(t, w));

    const double res = phib;
    const double objective = .5 * (res * res);

    const double opt_res = phib * alpha * fabs(c);

    sleqp_log_debug("Iteration %d, objective %e, residuum %e",
                    iteration,
                    objective,
                    opt_res);

    if(opt_res <= rel_tol)
    {
      break;
    }
  }

  sleqp_log_debug("LSQR solver terminated with an interior solution after %d iterations",
                  iteration);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_lsqr_solver_compute_step(SleqpLSQRSolver* solver,
                                             SleqpSparseVec* step)
{
  const double zero_eps = sleqp_params_get(solver->params,
                                           SLEQP_PARAM_ZERO_EPS);

  SleqpSparseVec* initial_step = sleqp_working_step_get_step(solver->working_step);

  SLEQP_CALL(lsqr_loop(solver));

  SLEQP_CALL(sleqp_sparse_vector_add(initial_step,
                                     solver->lsqr.x,
                                     zero_eps,
                                     step));

#ifndef NDEBUG

  if(sleqp_working_step_in_working_set(solver->working_step))
  {
    // Direction must be in working set
    bool in_working_set = false;

    const double eps = sleqp_params_get(solver->params,
                                        SLEQP_PARAM_EPS);

    SLEQP_NUM_ASSERT_PARAM(eps);

    SLEQP_CALL(sleqp_direction_in_working_set(solver->problem,
                                              solver->iterate,
                                              step,
                                              solver->dense_cache,
                                              eps,
                                              &in_working_set));

    sleqp_num_assert(in_working_set);
  }

#endif

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE lsqr_solver_free(SleqpLSQRSolver** star)
{
  SleqpLSQRSolver* solver = *star;

  if(!solver)
  {
    return SLEQP_OKAY;
  }

  {
    SLEQP_CALL(sleqp_sparse_vector_free(&solver->adjoint.product));

    SLEQP_CALL(sleqp_sparse_vector_free(&solver->adjoint.cons_product));
    SLEQP_CALL(sleqp_sparse_vector_free(&solver->adjoint.lsq_product));

    SLEQP_CALL(sleqp_sparse_vector_free(&solver->adjoint.cons_direction));
    SLEQP_CALL(sleqp_sparse_vector_free(&solver->adjoint.lsq_direction));
  }

  {
    SLEQP_CALL(sleqp_sparse_vector_free(&solver->forward.projected_direction));
    SLEQP_CALL(sleqp_sparse_vector_free(&solver->forward.lsq_product));
    SLEQP_CALL(sleqp_sparse_vector_free(&solver->forward.scaled_cons_product));
  }

  {
    SLEQP_CALL(sleqp_sparse_vector_free(&solver->lsqr.d));
    SLEQP_CALL(sleqp_sparse_vector_free(&solver->lsqr.x));

    SLEQP_CALL(sleqp_sparse_vector_free(&solver->lsqr.q));
    SLEQP_CALL(sleqp_sparse_vector_free(&solver->lsqr.p));

    SLEQP_CALL(sleqp_sparse_vector_free(&solver->lsqr.w));

    SLEQP_CALL(sleqp_sparse_vector_free(&solver->lsqr.v));
    SLEQP_CALL(sleqp_sparse_vector_free(&solver->lsqr.u));
  }

  SLEQP_CALL(sleqp_aug_jacobian_release(&solver->jacobian));

  sleqp_free(&solver->dense_cache);

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->sparse_cache));

  sleqp_free(&solver->removed_cons);

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->scaled_cons_residuals));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->linear_cons_residuals));

  SLEQP_CALL(sleqp_sparse_matrix_release(&solver->scaled_violated_cons_jac));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->linear_cons_val));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->rhs));
  SLEQP_CALL(sleqp_sparse_vector_free(&solver->lsq_rhs));
  SLEQP_CALL(sleqp_sparse_vector_free(&solver->lsq_residuals));

  SLEQP_CALL(sleqp_iterate_release(&solver->iterate));

  SLEQP_CALL(sleqp_params_release(&solver->params));
  SLEQP_CALL(sleqp_working_step_release(&solver->working_step));
  SLEQP_CALL(sleqp_problem_release(&solver->problem));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_lsqr_solver_capture(SleqpLSQRSolver* solver)
{
  ++solver->refcount;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_lsqr_solver_release(SleqpLSQRSolver** star)
{
  SleqpLSQRSolver* solver = *star;

  if(!solver)
  {
    return SLEQP_OKAY;
  }

  if(--solver->refcount == 0)
  {
    SLEQP_CALL(lsqr_solver_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
