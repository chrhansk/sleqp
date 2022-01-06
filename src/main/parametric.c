#include "parametric.h"

#include <math.h>

#include "cmp.h"
#include "fail.h"
#include "feas.h"
#include "log.h"
#include "mem.h"

struct SleqpParametricSolver
{
  int refcount;
  SleqpProblem* problem;
  SleqpParams* params;
  SleqpMerit* merit;
  SleqpLineSearch* linesearch;

  double exact_violation;

  double* cache;
  SleqpSparseVec* jacobian_product;
  SleqpSparseVec* combined_cons_val;
  SleqpSparseVec* last_hessian_direction;

  double trust_radius_increase;
  double trust_radius_decrease;
  int max_num_resolves;

  double penalty_parameter;
};

SLEQP_RETCODE
sleqp_parametric_solver_create(SleqpParametricSolver** star,
                               SleqpProblem* problem,
                               SleqpParams* params,
                               SleqpOptions* options,
                               SleqpMerit* merit,
                               SleqpLineSearch* linesearch)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpParametricSolver* solver = *star;

  solver->refcount = 1;

  solver->problem = problem;
  SLEQP_CALL(sleqp_problem_capture(solver->problem));

  solver->params = params;
  SLEQP_CALL(sleqp_params_capture(solver->params));

  solver->merit = merit;
  SLEQP_CALL(sleqp_merit_capture(solver->merit));

  solver->linesearch = linesearch;
  SLEQP_CALL(sleqp_linesearch_capture(solver->linesearch));

  const int num_variables   = sleqp_problem_num_vars(problem);
  const int num_constraints = sleqp_problem_num_cons(problem);

  SLEQP_CALL(sleqp_alloc_array(&solver->cache, num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->jacobian_product,
                                              num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->combined_cons_val,
                                              num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->last_hessian_direction,
                                              num_variables));

  SLEQP_PARAMETRIC_CAUCHY parametric_cauchy
    = sleqp_options_enum_value(options, SLEQP_OPTION_ENUM_PARAMETRIC_CAUCHY);

  if (parametric_cauchy == SLEQP_PARAMETRIC_CAUCHY_COARSE)
  {
    solver->trust_radius_increase = 2.;
    solver->trust_radius_decrease = .5;
    solver->max_num_resolves      = 5;
  }
  else
  {
    assert(parametric_cauchy == SLEQP_PARAMETRIC_CAUCHY_FINE);
    solver->trust_radius_increase = sqrt(2.);
    solver->trust_radius_decrease = sqrt(.5);
    solver->max_num_resolves      = 10;
  }

  return SLEQP_OKAY;
}

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_parametric_solver_set_penalty(SleqpParametricSolver* solver,
                                    double penalty_parameter)
{
  solver->penalty_parameter = penalty_parameter;
  return SLEQP_OKAY;
}

static SLEQP_RETCODE
has_sufficient_decrease(SleqpParametricSolver* solver,
                        SleqpIterate* iterate,
                        const SleqpSparseVec* direction,
                        SleqpSparseVec* hessian_direction,
                        const SleqpSparseVec* multipliers,
                        double* quadratic_merit,
                        bool* sufficient_decrease)
{
  double objective_dot, linear_violation, hessian_product;

  const double one = 1.;
  const double eta = sleqp_params_value(solver->params, SLEQP_PARAM_CAUCHY_ETA);

  SleqpProblem* problem = solver->problem;

  const int num_constraints = sleqp_problem_num_cons(problem);

  const double zero_eps
    = sleqp_params_value(solver->params, SLEQP_PARAM_ZERO_EPS);

  const double exact_violation = solver->exact_violation;

  const double penalty_parameter = solver->penalty_parameter;

  SLEQP_CALL(sleqp_sparse_vector_dot(sleqp_iterate_obj_grad(iterate),
                                     direction,
                                     &objective_dot));

  SLEQP_CALL(sleqp_sparse_matrix_vector_product(sleqp_iterate_cons_jac(iterate),
                                                direction,
                                                solver->cache));

  SLEQP_CALL(sleqp_sparse_vector_from_raw(solver->jacobian_product,
                                          solver->cache,
                                          num_constraints,
                                          zero_eps));

  SLEQP_CALL(sleqp_sparse_vector_add(sleqp_iterate_cons_val(iterate),
                                     solver->jacobian_product,
                                     zero_eps,
                                     solver->combined_cons_val));

  SLEQP_CALL(sleqp_violation_one_norm(problem,
                                      solver->combined_cons_val,
                                      &linear_violation));

  SLEQP_CALL(sleqp_problem_hess_prod(problem,
                                     &one,
                                     direction,
                                     multipliers,
                                     hessian_direction));

  SLEQP_CALL(
    sleqp_sparse_vector_dot(direction, hessian_direction, &hessian_product));

  *sufficient_decrease
    = ((penalty_parameter * (exact_violation - linear_violation)
        - objective_dot)
         * (1. - eta)
       >= (0.5 * hessian_product));

  *quadratic_merit = sleqp_iterate_obj_val(iterate) + objective_dot
                     + (penalty_parameter * linear_violation)
                     + .5 * hessian_product;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
search_forward(SleqpParametricSolver* solver,
               SleqpIterate* iterate,
               SleqpCauchy* cauchy_data,
               SleqpSparseVec* cauchy_direction,
               SleqpSparseVec* cauchy_hessian_direction,
               const SleqpSparseVec* multipliers,
               double* trust_radius,
               double* quadratic_merit)
{
  SleqpProblem* problem = solver->problem;

  const double eps = sleqp_params_value(solver->params, SLEQP_PARAM_EPS);

  const double penalty_parameter = solver->penalty_parameter;

  const double one = 1.;

  double last_quadratic_merit = *quadratic_merit;

  SLEQP_CALL(sleqp_sparse_vector_copy(cauchy_hessian_direction,
                                      solver->last_hessian_direction));

  for (int i = 0; i < solver->max_num_resolves; ++i)
  {
    sleqp_log_debug("Resolving with radius %.14e", *trust_radius);

    SLEQP_CALL(sleqp_cauchy_set_trust_radius(cauchy_data, (*trust_radius)));

    SLEQP_CALL(sleqp_cauchy_solve(cauchy_data,
                                  sleqp_iterate_obj_grad(iterate),
                                  penalty_parameter,
                                  SLEQP_CAUCHY_OBJECTIVE_TYPE_DEFAULT));

    SLEQP_CALL(sleqp_cauchy_get_direction(cauchy_data, cauchy_direction));

    {
      *quadratic_merit = 0.;

      SLEQP_CALL(sleqp_merit_linear(solver->merit,
                                    iterate,
                                    cauchy_direction,
                                    penalty_parameter,
                                    quadratic_merit));

      SLEQP_CALL(sleqp_problem_hess_prod(problem,
                                         &one,
                                         cauchy_direction,
                                         multipliers,
                                         cauchy_hessian_direction));

      double hessian_dot;

      SLEQP_CALL(sleqp_sparse_vector_dot(cauchy_direction,
                                         cauchy_hessian_direction,
                                         &hessian_dot));

      *quadratic_merit += .5 * hessian_dot;
    }

    if (sleqp_is_lt(*quadratic_merit, last_quadratic_merit, eps))
    {
      (*trust_radius) *= solver->trust_radius_increase;
    }
    else
    {
      // track back...
      (*trust_radius) /= solver->trust_radius_increase;

      *quadratic_merit = last_quadratic_merit;

      sleqp_log_debug(
        "Accepting trust radius of %.14e with a quadratic merit of %.14e",
        *trust_radius,
        last_quadratic_merit);

      SLEQP_CALL(sleqp_cauchy_set_trust_radius(cauchy_data, (*trust_radius)));

      SLEQP_CALL(sleqp_cauchy_solve(cauchy_data,
                                    sleqp_iterate_obj_grad(iterate),
                                    penalty_parameter,
                                    SLEQP_CAUCHY_OBJECTIVE_TYPE_DEFAULT));

      SLEQP_CALL(sleqp_cauchy_get_direction(cauchy_data, cauchy_direction));

      SLEQP_CALL(sleqp_sparse_vector_copy(solver->last_hessian_direction,
                                          cauchy_hessian_direction));

      break;
    }

    last_quadratic_merit = *quadratic_merit;

    SLEQP_CALL(sleqp_sparse_vector_copy(cauchy_hessian_direction,
                                        solver->last_hessian_direction));
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
search_backtracking(SleqpParametricSolver* solver,
                    SleqpIterate* iterate,
                    SleqpCauchy* cauchy_data,
                    SleqpSparseVec* cauchy_direction,
                    SleqpSparseVec* cauchy_hessian_direction,
                    const SleqpSparseVec* multipliers,
                    double* trust_radius,
                    double* quadratic_merit)
{
  double penalty_parameter = solver->penalty_parameter;

  int i = 0;

  for (; i < solver->max_num_resolves; ++i)
  {
    sleqp_log_debug("Resolving with radius %.14e", *trust_radius);

    SLEQP_CALL(sleqp_cauchy_set_trust_radius(cauchy_data, (*trust_radius)));

    SLEQP_CALL(sleqp_cauchy_solve(cauchy_data,
                                  sleqp_iterate_obj_grad(iterate),
                                  penalty_parameter,
                                  SLEQP_CAUCHY_OBJECTIVE_TYPE_DEFAULT));

    SLEQP_CALL(sleqp_cauchy_get_direction(cauchy_data, cauchy_direction));

#ifndef NDEBUG
    {
      const double eps = sleqp_params_value(solver->params, SLEQP_PARAM_EPS);

      const double step_norm = sleqp_sparse_vector_inf_norm(cauchy_direction);

      SLEQP_NUM_ASSERT_PARAM(eps);
      SLEQP_NUM_ASSERT_PARAM(step_norm);

      sleqp_num_assert(sleqp_is_leq(step_norm, *trust_radius, eps));
    }
#endif

    bool sufficient_decrease;

    SLEQP_CALL(has_sufficient_decrease(solver,
                                       iterate,
                                       cauchy_direction,
                                       cauchy_hessian_direction,
                                       multipliers,
                                       quadratic_merit,
                                       &sufficient_decrease));

    if (sufficient_decrease)
    {
      sleqp_log_debug(
        "Accepting radius %.14e, which provides sufficient decrease",
        *trust_radius);

      break;
    }
    else
    {
      (*trust_radius) *= solver->trust_radius_decrease;
    }
  }

  if (i == solver->max_num_resolves)
  {
    bool full_step;

    SLEQP_CALL(sleqp_linesearch_cauchy_step(solver->linesearch,
                                            cauchy_direction,
                                            multipliers,
                                            cauchy_hessian_direction,
                                            &full_step,
                                            quadratic_merit));
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_parametric_solver_solve(SleqpParametricSolver* solver,
                              SleqpIterate* iterate,
                              SleqpCauchy* cauchy_data,
                              SleqpSparseVec* cauchy_direction,
                              SleqpSparseVec* cauchy_hessian_direction,
                              const SleqpSparseVec* multipliers,
                              double* trust_radius,
                              double* quadratic_merit_value)
{
  SLEQP_CALL(sleqp_violation_one_norm(solver->problem,
                                      sleqp_iterate_cons_val(iterate),
                                      &solver->exact_violation));

  bool sufficient_decrease;

  SLEQP_CALL(has_sufficient_decrease(solver,
                                     iterate,
                                     cauchy_direction,
                                     cauchy_hessian_direction,
                                     multipliers,
                                     quadratic_merit_value,
                                     &sufficient_decrease));

  sleqp_log_debug("Beginning parametric solve, initial trust radius: %.14e, "
                  "quadratic merit: %.14e",
                  *trust_radius,
                  *quadratic_merit_value);

  if (sufficient_decrease)
  {
    sleqp_log_debug(
      "Initial radius provides sufficient decrease, searching forward");

    (*trust_radius) *= solver->trust_radius_increase;

    SLEQP_CALL(search_forward(solver,
                              iterate,
                              cauchy_data,
                              cauchy_direction,
                              cauchy_hessian_direction,
                              multipliers,
                              trust_radius,
                              quadratic_merit_value));
  }
  else
  {
    sleqp_log_debug(
      "Initial radius does not provide sufficient decrease, tracking back");

    (*trust_radius) *= solver->trust_radius_decrease;

    SLEQP_CALL(search_backtracking(solver,
                                   iterate,
                                   cauchy_data,
                                   cauchy_direction,
                                   cauchy_hessian_direction,
                                   multipliers,
                                   trust_radius,
                                   quadratic_merit_value));
  }

#ifndef NDEBUG
  {
    const double eps = sleqp_params_value(solver->params, SLEQP_PARAM_EPS);

    const double step_norm = sleqp_sparse_vector_inf_norm(cauchy_direction);

    SLEQP_NUM_ASSERT_PARAM(eps);
    SLEQP_NUM_ASSERT_PARAM(step_norm);

    sleqp_num_assert(sleqp_is_leq(step_norm, *trust_radius, eps));
  }
#endif

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_parametric_solver_capture(SleqpParametricSolver* solver)
{
  ++solver->refcount;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
parametric_solver_free(SleqpParametricSolver** star)
{
  SleqpParametricSolver* solver = *star;

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->last_hessian_direction));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->combined_cons_val));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->jacobian_product));

  SLEQP_CALL(sleqp_linesearch_release(&solver->linesearch));

  SLEQP_CALL(sleqp_merit_release(&solver->merit));

  SLEQP_CALL(sleqp_params_release(&solver->params));

  SLEQP_CALL(sleqp_problem_release(&solver->problem));

  sleqp_free(&solver->cache);

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_parametric_solver_release(SleqpParametricSolver** star)
{
  SleqpParametricSolver* solver = *star;

  if (!solver)
  {
    return SLEQP_OKAY;
  }

  if (--solver->refcount == 0)
  {
    SLEQP_CALL(parametric_solver_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
