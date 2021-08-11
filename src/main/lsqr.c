#include "lsqr.h"

#include <math.h>

#include "aug_jacobian.h"
#include "cmp.h"
#include "feas.h"
#include "lsq.h"
#include "mem.h"
#include "working_step.h"

#include "tr/tr_util.h"

typedef struct LSQR {
  SleqpSparseVec* u;
  SleqpSparseVec* v;
  SleqpSparseVec* w;

  SleqpSparseVec* p;
  SleqpSparseVec* q;

  SleqpSparseVec* x;

  SleqpSparseVec* d;

  SleqpSparseVec* x_proj;
  SleqpSparseVec* x_proj_prev;
} LSQR;

typedef struct Adjoint {
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
  SleqpSparseVec* lsq_step_forward;
  SleqpSparseVec* lsq_rhs;
  SleqpSparseVec* rhs;

  SleqpSparseVec* linear_cons_val;
  SleqpSparseVec* linear_cons_residuals;

  SleqpSparseMatrix* scaled_violated_cons_jac;
  SleqpSparseVec* scaled_cons_residuals;
  SleqpSparseVec* scaled_cons_product;

  int* removed_cons;
  int num_removed_cons;

  SleqpSparseVec* sparse_cache;
  double* dense_cache;

  SleqpSparseVec* projected_direction;

  SleqpAugJacobian* jacobian;

  double trust_radius;
  double penalty_parameter;

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

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->lsq_residuals,
                                              num_residuals));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->lsq_step_forward,
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

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->scaled_cons_product,
                                              0));

  SLEQP_CALL(sleqp_alloc_array(&solver->removed_cons,
                               num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->sparse_cache,
                                              num_constraints));

  SLEQP_CALL(sleqp_alloc_array(&solver->dense_cache,
                               num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->projected_direction,
                                              num_variables));

  {
    SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->lsqr.u,
                                                0));

    SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->lsqr.v,
                                                num_variables));

    SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->lsqr.w,
                                                num_variables));

    SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->lsqr.p,
                                                0));

    SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->lsqr.q,
                                                num_variables));

    SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->lsqr.x,
                                                num_variables));

    SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->lsqr.d,
                                                num_variables));

    SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->lsqr.x_proj,
                                                num_variables));

    SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->lsqr.x_proj_prev,
                                                num_variables));
  }

  {
    SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->adjoint.lsq_direction,
                                                num_residuals));

    SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->adjoint.cons_direction,
                                                0));

    SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->adjoint.lsq_product,
                                                num_variables));

    SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->adjoint.cons_product,
                                                num_variables));

    SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->adjoint.product,
                                                num_variables));
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
                                        solver->lsq_step_forward));

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(solver->lsq_residuals,
                                            solver->lsq_step_forward,
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
  SleqpFunc* func = sleqp_problem_func(problem);

  const double zero_eps = sleqp_params_get(solver->params,
                                           SLEQP_PARAM_ZERO_EPS);

  const int num_constraints = sleqp_problem_num_constraints(problem);
  const int num_variables = sleqp_problem_num_variables(problem);
  const int num_residuals = sleqp_lsq_func_num_residuals(func);

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

  const int num_required_cons = num_constraints - solver->num_removed_cons;

  SLEQP_CALL(sleqp_sparse_matrix_resize(solver->scaled_violated_cons_jac,
                                        num_required_cons,
                                        num_variables));

  SLEQP_CALL(sleqp_sparse_matrix_remove_rows(cons_jac,
                                             solver->scaled_violated_cons_jac,
                                             solver->removed_cons,
                                             solver->num_removed_cons));

  SLEQP_CALL(sleqp_sparse_matrix_scale(solver->scaled_violated_cons_jac,
                                       solver->penalty_parameter));

  SLEQP_CALL(sleqp_sparse_vector_clear(solver->scaled_cons_residuals));

  SLEQP_CALL(sleqp_sparse_vector_resize(solver->scaled_cons_residuals,
                                        num_required_cons));

  SLEQP_CALL(sleqp_sparse_vector_resize(solver->scaled_cons_product,
                                        num_required_cons));

  SLEQP_CALL(sleqp_sparse_vector_resize(solver->adjoint.cons_direction,
                                        num_required_cons));

  SLEQP_CALL(sleqp_sparse_vector_resize(solver->lsqr.u,
                                        num_residuals + num_required_cons));

  SLEQP_CALL(sleqp_sparse_vector_resize(solver->lsqr.p,
                                        num_residuals + num_required_cons));

  SLEQP_CALL(sleqp_sparse_vector_remove_entries(solver->linear_cons_residuals,
                                                solver->scaled_cons_residuals,
                                                solver->removed_cons,
                                                solver->num_removed_cons));

  SLEQP_CALL(sleqp_sparse_vector_scale(solver->scaled_cons_residuals,
                                       solver->penalty_parameter));

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE create_removed_cons(SleqpLSQRSolver* solver)
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
SLEQP_RETCODE compute_rhs(SleqpLSQRSolver* solver)
{
  SleqpProblem* problem = solver->problem;
  SleqpFunc* func = sleqp_problem_func(problem);

  const int num_residuals = sleqp_lsq_func_num_residuals(func);
  const int num_constraints = sleqp_problem_num_constraints(problem);
  const int num_required_cons = num_constraints - solver->num_removed_cons;

  SLEQP_CALL(compute_lsq_rhs(solver));

  SLEQP_CALL(compute_cons_rhs(solver));

  SLEQP_CALL(sleqp_sparse_vector_clear(solver->rhs));

  SLEQP_CALL(sleqp_sparse_vector_resize(solver->rhs,
                                        num_residuals + num_required_cons));

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

  SLEQP_CALL(create_removed_cons(solver));

  solver->trust_radius = sleqp_working_step_get_reduced_trust_radius(solver->working_step);
  solver->penalty_parameter = penalty_parameter;

  SLEQP_CALL(compute_rhs(solver));

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE forward_product(SleqpLSQRSolver* solver,
                              SleqpSparseVec* direction,
                              SleqpSparseVec* product)
{
  SleqpProblem* problem = solver->problem;
  SleqpFunc* func = sleqp_problem_func(problem);

  const double zero_eps = sleqp_params_get(solver->params,
                                           SLEQP_PARAM_ZERO_EPS);

  const int num_constraints = sleqp_problem_num_constraints(problem);
  const int num_required_cons = num_constraints - solver->num_removed_cons;

  SLEQP_CALL(sleqp_aug_jacobian_projection(solver->jacobian,
                                           direction,
                                           solver->projected_direction,
                                           NULL));

  SLEQP_CALL(sleqp_lsq_func_jac_forward(func,
                                        direction,
                                        solver->lsq_step_forward));

  SLEQP_CALL(sleqp_sparse_matrix_vector_product(solver->scaled_violated_cons_jac,
                                                solver->projected_direction,
                                                solver->dense_cache));

  SLEQP_CALL(sleqp_sparse_vector_from_raw(solver->scaled_cons_product,
                                          solver->dense_cache,
                                          num_required_cons,
                                          zero_eps));

  SLEQP_CALL(sleqp_sparse_vector_concat(solver->lsq_step_forward,
                                        solver->scaled_cons_product,
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
  const int num_constraints = sleqp_problem_num_constraints(problem);
  const int num_required_cons = num_constraints - solver->num_removed_cons;

  // split product
  {
    SLEQP_CALL(sleqp_sparse_vector_clear(solver->adjoint.lsq_direction));
    SLEQP_CALL(sleqp_sparse_vector_reserve(solver->adjoint.lsq_direction,
                                           SLEQP_MIN(num_residuals, direction->nnz)));

    SLEQP_CALL(sleqp_sparse_vector_clear(solver->adjoint.cons_direction));
    SLEQP_CALL(sleqp_sparse_vector_reserve(solver->adjoint.cons_direction,
                                           SLEQP_MIN(num_required_cons, direction->nnz)));

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
    SLEQP_CALL(sleqp_sparse_vector_scale(vec,
                                         1./(*norm)));
  }

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE lsqr_loop(SleqpLSQRSolver* solver,
                        SleqpSparseVec* step)
{
  const double eps = sleqp_params_get(solver->params,
                                      SLEQP_PARAM_EPS);

  const double zero_eps = sleqp_params_get(solver->params,
                                           SLEQP_PARAM_ZERO_EPS);

  const int n = sleqp_problem_num_variables(solver->problem);
  SleqpSparseVec* x = solver->lsqr.x;
  SleqpSparseVec* u = solver->lsqr.u;
  SleqpSparseVec* v = solver->lsqr.v;
  SleqpSparseVec* w = solver->lsqr.w;

  SleqpSparseVec* p = solver->lsqr.p;
  SleqpSparseVec* q = solver->lsqr.q;

  SleqpSparseVec* d = solver->lsqr.d;

  SleqpSparseVec* b = solver->rhs;
  SleqpSparseVec* t = solver->sparse_cache;

  SleqpSparseVec* x_proj = solver->lsqr.x_proj;
  SleqpSparseVec* x_proj_prev = solver->lsqr.x_proj_prev;

  double alpha, beta;

  SLEQP_CALL(sleqp_sparse_vector_copy(b, u));
  SLEQP_CALL(normalize(u, &beta));

  SLEQP_CALL(adjoint_product(solver, u, v));
  SLEQP_CALL(normalize(v, &alpha));

  SLEQP_CALL(sleqp_sparse_vector_copy(v, w));

  SLEQP_CALL(sleqp_sparse_vector_clear(x));
  SLEQP_CALL(sleqp_sparse_vector_clear(x_proj));

  double phib = beta;
  double rhob = alpha;

  for(int i = 0; i < n; ++i)
  {
    SLEQP_CALL(forward_product(solver, v, p));
    SLEQP_CALL(sleqp_sparse_vector_add_scaled(p, u, 1., -alpha, zero_eps, t));
    SLEQP_CALL(sleqp_sparse_vector_copy(t, u));
    SLEQP_CALL(normalize(u, &beta));

    SLEQP_CALL(adjoint_product(solver, u, q));
    SLEQP_CALL(sleqp_sparse_vector_add_scaled(q, v, 1., -beta, zero_eps, t));
    SLEQP_CALL(sleqp_sparse_vector_copy(t, v));
    SLEQP_CALL(normalize(v, &alpha));

    double rho = hypot(rhob, beta);
    const double c = rhob / rho;
    const double s = beta / rho;
    double theta = s * alpha;
    rhob = -c*alpha;
    double phi = c*phib;
    phib = s*phib;

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(x,
                                              w,
                                              1.,
                                              phi / rho,
                                              zero_eps,
                                              t));

    SLEQP_CALL(sleqp_sparse_vector_copy(x_proj, x_proj_prev));

    SLEQP_CALL(sleqp_sparse_vector_copy(t, x));

    SLEQP_CALL(sleqp_aug_jacobian_projection(solver->jacobian,
                                             x,
                                             x_proj,
                                             NULL));

    const double norm = sleqp_sparse_vector_norm(x_proj);

    if(sleqp_is_gt(norm, solver->trust_radius, eps))
    {
      SLEQP_CALL(sleqp_sparse_vector_add_scaled(x_proj_prev,
                                                x_proj,
                                                -1.,
                                                1.,
                                                zero_eps,
                                                d));

      SLEQP_CALL(sleqp_tr_compute_bdry_sol(x_proj_prev,
                                           d,
                                           solver->params,
                                           solver->trust_radius,
                                           step));


      return SLEQP_OKAY;
    }

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(v,
                                              w,
                                              1.,
                                              - theta / rho,
                                              zero_eps,
                                              t));

    SLEQP_CALL(sleqp_sparse_vector_copy(t, w));
  }

  SLEQP_CALL(sleqp_sparse_vector_copy(x_proj, step));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_lsqr_solver_compute_step(SleqpLSQRSolver* solver,
                                             SleqpSparseVec* lsqr_step)
{
  SLEQP_CALL(lsqr_loop(solver, lsqr_step));

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
    SLEQP_CALL(sleqp_sparse_vector_free(&solver->lsqr.x_proj_prev));
    SLEQP_CALL(sleqp_sparse_vector_free(&solver->lsqr.x_proj));
    SLEQP_CALL(sleqp_sparse_vector_free(&solver->lsqr.d));
    SLEQP_CALL(sleqp_sparse_vector_free(&solver->lsqr.x));

    SLEQP_CALL(sleqp_sparse_vector_free(&solver->lsqr.q));
    SLEQP_CALL(sleqp_sparse_vector_free(&solver->lsqr.p));

    SLEQP_CALL(sleqp_sparse_vector_free(&solver->lsqr.w));

    SLEQP_CALL(sleqp_sparse_vector_free(&solver->lsqr.v));
    SLEQP_CALL(sleqp_sparse_vector_free(&solver->lsqr.u));
  }

  SLEQP_CALL(sleqp_aug_jacobian_release(&solver->jacobian));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->projected_direction));

  sleqp_free(&solver->dense_cache);

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->sparse_cache));

  sleqp_free(&solver->removed_cons);

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->scaled_cons_product));
  SLEQP_CALL(sleqp_sparse_vector_free(&solver->scaled_cons_residuals));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->linear_cons_residuals));

  SLEQP_CALL(sleqp_sparse_matrix_release(&solver->scaled_violated_cons_jac));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->linear_cons_val));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->rhs));
  SLEQP_CALL(sleqp_sparse_vector_free(&solver->lsq_rhs));
  SLEQP_CALL(sleqp_sparse_vector_free(&solver->lsq_step_forward));
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
