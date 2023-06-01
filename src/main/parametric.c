#include "parametric.h"

#include <math.h>

#include "cmp.h"
#include "direction.h"
#include "fail.h"
#include "feas.h"
#include "log.h"
#include "mem.h"
#include "sparse/pub_vec.h"

struct SleqpParametricSolver
{
  int refcount;
  SleqpProblem* problem;
  SleqpParams* params;
  SleqpMerit* merit;
  SleqpLineSearch* linesearch;

  double exact_violation;

  double* cache;
  SleqpVec* jacobian_product;
  SleqpVec* combined_cons_val;

  SleqpDirection* last_direction;

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

  const int num_constraints = sleqp_problem_num_cons(problem);

  SLEQP_CALL(sleqp_alloc_array(&solver->cache, num_constraints));

  SLEQP_CALL(
    sleqp_vec_create_empty(&solver->jacobian_product, num_constraints));

  SLEQP_CALL(
    sleqp_vec_create_empty(&solver->combined_cons_val, num_constraints));

  SLEQP_CALL(sleqp_direction_create(&solver->last_direction, problem, params));

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
                        SleqpDirection* direction,
                        double* quadratic_merit,
                        bool* sufficient_decrease)
{
  double linear_violation, hessian_product;

  const double eta = sleqp_params_value(solver->params, SLEQP_PARAM_CAUCHY_ETA);

  SleqpProblem* problem = solver->problem;

  const double zero_eps
    = sleqp_params_value(solver->params, SLEQP_PARAM_ZERO_EPS);

  const double exact_violation = solver->exact_violation;

  const double penalty_parameter = solver->penalty_parameter;

  const double objective_dot = *sleqp_direction_obj_grad(direction);

  SleqpVec* direction_cons_jac = sleqp_direction_cons_jac(direction);

  SLEQP_CALL(sleqp_vec_add(sleqp_iterate_cons_val(iterate),
                           direction_cons_jac,
                           zero_eps,
                           solver->combined_cons_val));

  SLEQP_CALL(sleqp_total_violation(problem,
                                   solver->combined_cons_val,
                                   &linear_violation));

  SLEQP_CALL(sleqp_vec_dot(sleqp_direction_primal(direction),
                           sleqp_direction_hess(direction),
                           &hessian_product));

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
               SleqpDirection* cauchy_direction,
               const SleqpVec* multipliers,
               double* trust_radius,
               double* quadratic_merit)
{
  SleqpProblem* problem = solver->problem;

  const double eps = sleqp_params_value(solver->params, SLEQP_PARAM_EPS);

  const double penalty_parameter = solver->penalty_parameter;

  double last_quadratic_merit = *quadratic_merit;

  SleqpVec* direction_primal = sleqp_direction_primal(cauchy_direction);
  SleqpVec* direction_hess   = sleqp_direction_hess(cauchy_direction);

  SLEQP_CALL(sleqp_direction_copy(cauchy_direction, solver->last_direction));

  for (int i = 0; i < solver->max_num_resolves; ++i)
  {
    sleqp_log_debug("Resolving with radius %.14e", *trust_radius);

    SLEQP_CALL(sleqp_cauchy_set_trust_radius(cauchy_data, (*trust_radius)));

    SLEQP_CALL(sleqp_cauchy_solve(cauchy_data,
                                  sleqp_iterate_obj_grad(iterate),
                                  penalty_parameter,
                                  SLEQP_CAUCHY_OBJTYPE_DEFAULT));

    SLEQP_CALL(sleqp_cauchy_lp_step(cauchy_data, direction_primal));

    {
      *quadratic_merit = 0.;

      SLEQP_CALL(sleqp_merit_linear(solver->merit,
                                    iterate,
                                    cauchy_direction,
                                    penalty_parameter,
                                    quadratic_merit));

      SLEQP_CALL(sleqp_problem_hess_prod(problem,
                                         direction_primal,
                                         multipliers,
                                         direction_hess));

      double hessian_dot;

      SLEQP_CALL(sleqp_vec_dot(direction_primal, direction_hess, &hessian_dot));

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

      // TODO: Do we need to resolve?
      SLEQP_CALL(sleqp_cauchy_solve(cauchy_data,
                                    sleqp_iterate_obj_grad(iterate),
                                    penalty_parameter,
                                    SLEQP_CAUCHY_OBJTYPE_DEFAULT));

      SLEQP_CALL(
        sleqp_direction_copy(solver->last_direction, cauchy_direction));

      break;
    }

    last_quadratic_merit = *quadratic_merit;

    SLEQP_CALL(sleqp_direction_copy(cauchy_direction, solver->last_direction));
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
search_backtracking(SleqpParametricSolver* solver,
                    SleqpIterate* iterate,
                    SleqpCauchy* cauchy_data,
                    SleqpDirection* cauchy_direction,
                    const SleqpVec* multipliers,
                    double* trust_radius,
                    double* quadratic_merit)
{
  double penalty_parameter = solver->penalty_parameter;

  int i = 0;

  SleqpVec* direction_primal = sleqp_direction_primal(cauchy_direction);

  for (; i < solver->max_num_resolves; ++i)
  {
    sleqp_log_debug("Resolving with radius %.14e", *trust_radius);

    SLEQP_CALL(sleqp_cauchy_set_trust_radius(cauchy_data, (*trust_radius)));

    SLEQP_CALL(sleqp_cauchy_solve(cauchy_data,
                                  sleqp_iterate_obj_grad(iterate),
                                  penalty_parameter,
                                  SLEQP_CAUCHY_OBJTYPE_DEFAULT));

    SLEQP_CALL(sleqp_cauchy_lp_step(cauchy_data, direction_primal));

    const double zero_eps
      = sleqp_params_value(solver->params, SLEQP_PARAM_ZERO_EPS);

    SLEQP_CALL(sleqp_direction_reset(cauchy_direction,
                                     solver->problem,
                                     iterate,
                                     multipliers,
                                     solver->cache,
                                     zero_eps));

#if SLEQP_DEBUG
    {
      const double eps = sleqp_params_value(solver->params, SLEQP_PARAM_EPS);

      const double step_norm = sleqp_vec_inf_norm(direction_primal);

      SLEQP_NUM_ASSERT_PARAM(eps);
      SLEQP_NUM_ASSERT_PARAM(step_norm);

      sleqp_num_assert(sleqp_is_leq(step_norm, *trust_radius, eps));
    }
#endif

    bool sufficient_decrease;

    SLEQP_CALL(has_sufficient_decrease(solver,
                                       iterate,
                                       cauchy_direction,
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
                                            &full_step,
                                            quadratic_merit));
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_parametric_solver_solve(SleqpParametricSolver* solver,
                              SleqpIterate* iterate,
                              SleqpCauchy* cauchy_data,
                              const SleqpVec* lp_step,
                              const SleqpVec* multipliers,
                              SleqpDirection* cauchy_direction,
                              double* trust_radius,
                              double* quadratic_merit_value)
{
  {
    SleqpVec* direction_primal = sleqp_direction_primal(cauchy_direction);
    SLEQP_CALL(sleqp_vec_copy(lp_step, direction_primal));

    const double zero_eps
      = sleqp_params_value(solver->params, SLEQP_PARAM_ZERO_EPS);

    SLEQP_CALL(sleqp_direction_reset(cauchy_direction,
                                     solver->problem,
                                     iterate,
                                     multipliers,
                                     solver->cache,
                                     zero_eps));
  }

  SLEQP_CALL(sleqp_total_violation(solver->problem,
                                   sleqp_iterate_cons_val(iterate),
                                   &solver->exact_violation));

  bool sufficient_decrease;

  SLEQP_CALL(has_sufficient_decrease(solver,
                                     iterate,
                                     cauchy_direction,
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
                                   multipliers,
                                   trust_radius,
                                   quadratic_merit_value));
  }

#if SLEQP_DEBUG
  {
    const double eps = sleqp_params_value(solver->params, SLEQP_PARAM_EPS);

    SleqpVec* cauchy_step = sleqp_direction_primal(cauchy_direction);

    const double step_norm = sleqp_vec_inf_norm(cauchy_step);

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

  SLEQP_CALL(sleqp_direction_release(&solver->last_direction));

  SLEQP_CALL(sleqp_vec_free(&solver->combined_cons_val));

  SLEQP_CALL(sleqp_vec_free(&solver->jacobian_product));

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
