#include "linesearch.h"

#include "cmp.h"
#include "direction.h"
#include "fail.h"
#include "feas.h"
#include "mem.h"
#include "pub_settings.h"
#include "sparse/pub_vec.h"
#include "util.h"

#define LINESEARCH_MAX_IT 10000

typedef struct BreakPoint
{
  double slope_change;
  double point;
} BreakPoint;

static int
compare_breakpoints(const void* first, const void* second)
{
  const double first_point  = ((const BreakPoint*)first)->point;
  const double second_point = ((const BreakPoint*)second)->point;

  if (first_point < second_point)
  {
    return -1;
  }
  else if (first_point > second_point)
  {
    return 1;
  }

  return 0;
}

struct SleqpLineSearch
{
  int refcount;

  SleqpProblem* problem;
  SleqpSettings* settings;
  SleqpMerit* merit;

  SleqpIterate* iterate;
  double penalty_parameter;
  double trust_radius;

  SleqpVec* cauchy_point;

  SleqpVec* cauchy_newton_jacobian_prod;
  SleqpVec* cauchy_cons_val;

  SleqpVec* cauchy_violation;

  SleqpVec* combined_cons_val;

  SleqpVec* cauchy_newton_direction;

  SleqpVec* violated_multipliers;

  SleqpDirection* test_direction;
  double* cache;

  int num_breakpoints;
  int max_num_breakpoints;
  BreakPoint* breakpoints;

  SleqpTimer* timer;
};

SLEQP_RETCODE
sleqp_linesearch_create(SleqpLineSearch** star,
                        SleqpProblem* problem,
                        SleqpSettings* settings,
                        SleqpMerit* merit)
{
  SLEQP_CALL(sleqp_malloc(star));

  const int num_variables = sleqp_problem_num_vars(problem);

  const int num_constraints = sleqp_problem_num_cons(problem);

  SleqpLineSearch* linesearch = *star;

  *linesearch = (SleqpLineSearch){0};

  linesearch->refcount = 1;

  linesearch->problem = problem;
  SLEQP_CALL(sleqp_problem_capture(linesearch->problem));

  SLEQP_CALL(sleqp_settings_capture(settings));
  linesearch->settings = settings;

  SLEQP_CALL(sleqp_merit_capture(merit));

  linesearch->merit = merit;

  SLEQP_CALL(sleqp_vec_create_empty(&linesearch->cauchy_point, num_variables));

  SLEQP_CALL(sleqp_vec_create_empty(&linesearch->cauchy_newton_jacobian_prod,
                                    num_constraints));

  SLEQP_CALL(
    sleqp_vec_create_empty(&linesearch->cauchy_cons_val, num_constraints));

  SLEQP_CALL(
    sleqp_vec_create_empty(&linesearch->cauchy_violation, num_constraints));

  SLEQP_CALL(
    sleqp_vec_create_empty(&linesearch->combined_cons_val, num_constraints));

  SLEQP_CALL(sleqp_vec_create_empty(&linesearch->cauchy_newton_direction,
                                    num_variables));

  SLEQP_CALL(
    sleqp_vec_create_empty(&linesearch->violated_multipliers, num_constraints));

  SLEQP_CALL(
    sleqp_direction_create(&linesearch->test_direction, problem, settings));

  SLEQP_CALL(sleqp_alloc_array(&linesearch->cache,
                               SLEQP_MAX(num_variables, num_constraints)));

  SLEQP_CALL(sleqp_timer_create(&linesearch->timer));

  linesearch->max_num_breakpoints = 0;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_linesearch_set_iterate(SleqpLineSearch* linesearch,
                             SleqpIterate* iterate,
                             double penalty_parameter,
                             double trust_radius)
{
  SLEQP_CALL(sleqp_iterate_release(&linesearch->iterate));

  linesearch->iterate = iterate;

  SLEQP_CALL(sleqp_iterate_capture(linesearch->iterate));

  linesearch->penalty_parameter = penalty_parameter;
  linesearch->trust_radius      = trust_radius;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_linesearch_cauchy_step(SleqpLineSearch* linesearch,
                             SleqpDirection* direction,
                             bool* full_step,
                             double* quadratic_merit_value)
{
  SleqpProblem* problem = linesearch->problem;

  SleqpIterate* iterate          = linesearch->iterate;
  const double penalty_parameter = linesearch->penalty_parameter;
  const double trust_radius      = linesearch->trust_radius;

  const double eps = sleqp_settings_real_value(linesearch->settings, SLEQP_SETTINGS_REAL_EPS);

  const double zero_eps
    = sleqp_settings_real_value(linesearch->settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  SLEQP_CALL(sleqp_timer_start(linesearch->timer));

  double exact_violation;

  SleqpVec* direction_primal   = sleqp_direction_primal(direction);
  SleqpVec* direction_cons_jac = sleqp_direction_cons_jac(direction);
  SleqpVec* direction_hess     = sleqp_direction_hess(direction);

  SLEQP_CALL(sleqp_total_violation(problem,
                                   sleqp_iterate_cons_val(iterate),
                                   &exact_violation));

  (*full_step)             = true;
  (*quadratic_merit_value) = 0.;

  const double exact_merit_value
    = sleqp_iterate_obj_val(iterate) + penalty_parameter * exact_violation;

  double hessian_product;

  SLEQP_CALL(sleqp_vec_dot(direction_primal, direction_hess, &hessian_product));

  double objective_dot;

  objective_dot = *sleqp_direction_obj_grad(direction);

  double delta = 1.;

  {
    double direction_norm = sleqp_vec_norm(direction_primal);

    double direction_factor = trust_radius / direction_norm;

    if (direction_factor < 1)
    {
      delta = direction_factor;

      hessian_product *= (delta * delta);
      objective_dot *= delta;

      SLEQP_CALL(sleqp_vec_scale(direction_cons_jac, delta));
    }
  }

  const double eta
    = sleqp_settings_real_value(linesearch->settings, SLEQP_SETTINGS_REAL_CAUCHY_ETA);
  const double tau
    = sleqp_settings_real_value(linesearch->settings, SLEQP_SETTINGS_REAL_CAUCHY_TAU);

  double linear_violation;

  int iteration = 0;

  sleqp_log_debug("Beginning Cauchy line search, exact merit value: %f",
                  exact_merit_value);

  bool zero_direction = false;

  for (iteration = 0; iteration < LINESEARCH_MAX_IT; ++iteration)
  {
    double linear_merit_value = 0.;

    // compute linear merit value
    {
      linear_merit_value += sleqp_iterate_obj_val(iterate) + objective_dot;

      SLEQP_CALL(sleqp_vec_add(sleqp_iterate_cons_val(iterate),
                               direction_cons_jac,
                               zero_eps,
                               linesearch->combined_cons_val));

      SLEQP_CALL(sleqp_total_violation(problem,
                                       linesearch->combined_cons_val,
                                       &linear_violation));

      linear_merit_value += penalty_parameter * linear_violation;
    }

    // compute quadratic merit value
    {
      (*quadratic_merit_value) = linear_merit_value + (0.5 * hessian_product);
    }

    sleqp_log_debug("Cauchy line search iteration %d, step length: %g, linear "
                    "merit value: %f, quadratic merit value: %f",
                    iteration,
                    delta,
                    linear_merit_value,
                    *quadratic_merit_value);

    // check condition
    // This is a numerically more sensible version of the condition
    // ((exact_merit_value - (*quadratic_merit_value)) >= eta*(exact_merit_value
    // - linear_merit_value))
    if ((penalty_parameter * (exact_violation - linear_violation)
         - objective_dot)
          * (1. - eta)
        >= (0.5 * hessian_product))
    {
      break;
    }

    // update products
    hessian_product *= (tau * tau);

    objective_dot *= tau;

    SLEQP_CALL(sleqp_vec_scale(direction_cons_jac, tau));

    delta *= tau;

    (*full_step) = false;

    if (sleqp_is_zero(delta, eps))
    {
      zero_direction = true;
      delta          = 0.;
      break;
    }
  }

  assert(iteration != LINESEARCH_MAX_IT);

  if (zero_direction)
  {
    SLEQP_CALL(sleqp_direction_set_zero(direction));
  }
  else
  {
    SLEQP_CALL(sleqp_vec_scale(direction_primal, delta));
    SLEQP_CALL(sleqp_vec_scale(direction_hess, delta));

    *sleqp_direction_obj_grad(direction) = objective_dot;
  }

  sleqp_log_debug("Cauchy line search terminated after %d iterations (step "
                  "length: %f, quadratic merit: %f)",
                  iteration,
                  delta,
                  (*quadratic_merit_value));

  sleqp_assert_is_leq(sleqp_vec_norm(direction_primal), trust_radius, eps);

  SLEQP_CALL(sleqp_timer_stop(linesearch->timer));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_linesearch_trial_step(SleqpLineSearch* linesearch,
                            const SleqpDirection* cauchy_direction,
                            const double cauchy_quadratic_merit_value,
                            const SleqpDirection* newton_direction,
                            const SleqpVec* multipliers,
                            SleqpDirection* trial_direction,
                            double* step_length,
                            double* trial_quadratic_merit_value)
{
  SleqpProblem* problem          = linesearch->problem;
  SleqpMerit* merit_data         = linesearch->merit;
  SleqpIterate* iterate          = linesearch->iterate;
  const double penalty_parameter = linesearch->penalty_parameter;

  const double eps = sleqp_settings_real_value(linesearch->settings, SLEQP_SETTINGS_REAL_EPS);

  SLEQP_NUM_ASSERT_PARAM(eps);

  const double zero_eps
    = sleqp_settings_real_value(linesearch->settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  SleqpVec* cauchy_step = sleqp_direction_primal(cauchy_direction);
  SleqpVec* newton_step = sleqp_direction_primal(newton_direction);

  SleqpVec* cauchy_hessian_step = sleqp_direction_hess(cauchy_direction);
  SleqpVec* newton_hessian_step = sleqp_direction_hess(newton_direction);

  SleqpVec* cauchy_cons_jac = sleqp_direction_cons_jac(cauchy_direction);
  SleqpVec* newton_cons_jac = sleqp_direction_cons_jac(newton_direction);

  SLEQP_CALL(sleqp_timer_start(linesearch->timer));

  // Compute Cauchy-Newton direction
  {
    SLEQP_CALL(sleqp_vec_add_scaled(newton_step,
                                    cauchy_step,
                                    1.,
                                    -1.,
                                    zero_eps,
                                    linesearch->cauchy_newton_direction));
  }

  // Compute initial scalar products for the linear merit function
  double cauchy_gradient_dot = *sleqp_direction_obj_grad(cauchy_direction);

  double newton_gradient_dot = *sleqp_direction_obj_grad(newton_direction);

  // Compute initial scalar products for the quadratic merit function
  double cauchy_cauchy_product;

  SLEQP_CALL(
    sleqp_vec_dot(cauchy_step, cauchy_hessian_step, &cauchy_cauchy_product));

  double cauchy_newton_product;

  SLEQP_CALL(
    sleqp_vec_dot(cauchy_step, newton_hessian_step, &cauchy_newton_product));

  double newton_newton_product;

  SLEQP_CALL(
    sleqp_vec_dot(newton_step, newton_hessian_step, &newton_newton_product));

  const double eta
    = sleqp_settings_real_value(linesearch->settings, SLEQP_SETTINGS_REAL_LINESEARCH_ETA);

  const double tau
    = sleqp_settings_real_value(linesearch->settings, SLEQP_SETTINGS_REAL_LINESEARCH_TAU);

  const double cutoff_threshold
    = sleqp_settings_real_value(linesearch->settings, SLEQP_SETTINGS_REAL_LINESEARCH_CUTOFF);

  int iteration = 0;

  double alpha = 1.;

  // compute the maximum trial step length
  {
    SLEQP_CALL(sleqp_vec_add(sleqp_iterate_primal(iterate),
                             cauchy_step,
                             zero_eps,
                             linesearch->cauchy_point));

    SLEQP_CALL(sleqp_max_step_length(linesearch->cauchy_point,
                                     linesearch->cauchy_newton_direction,
                                     sleqp_problem_vars_lb(problem),
                                     sleqp_problem_vars_ub(problem),
                                     &alpha));

    assert(alpha >= 0);
    assert(alpha <= 1);
  }

  double quadratic_merit_gradient_product = 0.;

  // compute merit gradient products
  {
    SLEQP_CALL(sleqp_vec_add(sleqp_iterate_cons_val(iterate),
                             cauchy_cons_jac,
                             zero_eps,
                             linesearch->cauchy_cons_val));

    // use the violated multipliers in 0, +/-1 to filter the Jacobian product
    SLEQP_CALL(sleqp_violated_cons_multipliers(problem,
                                               linesearch->cauchy_cons_val,
                                               linesearch->violated_multipliers,
                                               NULL));

    double jacobian_dot;

    SLEQP_CALL(sleqp_vec_dot(linesearch->violated_multipliers,
                             cauchy_cons_jac,
                             &jacobian_dot));

    double quadratic_merit_gradient_cauchy
      = cauchy_gradient_dot + jacobian_dot + cauchy_cauchy_product;

    SLEQP_CALL(sleqp_vec_dot(linesearch->violated_multipliers,
                             newton_cons_jac,
                             &jacobian_dot));

    double quadratic_merit_gradient_newton
      = newton_gradient_dot + jacobian_dot + cauchy_newton_product;

    quadratic_merit_gradient_product
      = (quadratic_merit_gradient_newton - quadratic_merit_gradient_cauchy);
  }

  if (alpha <= cutoff_threshold)
  {
    (*step_length)                 = 0.;
    (*trial_quadratic_merit_value) = cauchy_quadratic_merit_value;

    SLEQP_CALL(sleqp_direction_copy(cauchy_direction, trial_direction));

    SLEQP_CALL(sleqp_timer_stop(linesearch->timer));

    return SLEQP_OKAY;
  }

  sleqp_log_debug("Beginning Cauchy-Newton line search, initial step length: "
                  "%g, quadratic merit at Cauchy point: %g",
                  alpha,
                  cauchy_quadratic_merit_value);

  for (iteration = 0; iteration < LINESEARCH_MAX_IT; ++iteration)
  {
    double linear_merit_value    = 0.;
    double quadratic_merit_value = 0.;

    // Compute linear merit
    {
      linear_merit_value += sleqp_iterate_obj_val(iterate)
                            + ((1. - alpha) * cauchy_gradient_dot)
                            + ((alpha)*newton_gradient_dot);

      SLEQP_CALL(sleqp_vec_add_scaled(sleqp_iterate_cons_val(iterate),
                                      cauchy_cons_jac,
                                      1.,
                                      (1. - alpha),
                                      zero_eps,
                                      linesearch->cauchy_cons_val));

      SLEQP_CALL(sleqp_vec_add_scaled(linesearch->cauchy_cons_val,
                                      newton_cons_jac,
                                      1.,
                                      alpha,
                                      zero_eps,
                                      linesearch->combined_cons_val));

      double total_violation;

      SLEQP_CALL(sleqp_total_violation(problem,
                                       linesearch->combined_cons_val,
                                       &total_violation));

      linear_merit_value += penalty_parameter * total_violation;
    }

    // Compute quadratic term

    {
      double quadratic_term = 0.;

      quadratic_term
        += 0.5 * (1. - alpha) * (1. - alpha) * cauchy_cauchy_product;
      quadratic_term += alpha
                        * ((1. - alpha) * cauchy_newton_product
                           + 0.5 * alpha * newton_newton_product);

      quadratic_merit_value = linear_merit_value + quadratic_term;
    }

#if SLEQP_DEBUG

    {
      SleqpVec* test_step = sleqp_direction_primal(linesearch->test_direction);

      SLEQP_CALL(sleqp_vec_add_scaled(cauchy_step,
                                      newton_step,
                                      1. - alpha,
                                      alpha,
                                      zero_eps,
                                      test_step));

      SLEQP_CALL(sleqp_direction_reset(linesearch->test_direction,
                                       problem,
                                       iterate,
                                       multipliers,
                                       linesearch->cache,
                                       zero_eps));
    }

    {
      double actual_linear_merit;

      SLEQP_CALL(sleqp_merit_linear(merit_data,
                                    iterate,
                                    linesearch->test_direction,
                                    penalty_parameter,
                                    &actual_linear_merit));

      sleqp_assert_is_eq(linear_merit_value, actual_linear_merit, eps);
    }

    // Check quadratic merit
    {
      double actual_quadratic_merit;

      SLEQP_CALL(sleqp_merit_quadratic(merit_data,
                                       iterate,
                                       linesearch->test_direction,
                                       penalty_parameter,
                                       &actual_quadratic_merit));

      sleqp_assert_is_eq(quadratic_merit_value, actual_quadratic_merit, eps);
    }

#endif

    const double scaled_product = alpha * quadratic_merit_gradient_product;

    sleqp_log_debug("Cauchy-Newton line search iteration %d, step length: %g, "
                    "quadratic merit value: %g, scaled inner product: %g",
                    iteration,
                    alpha,
                    quadratic_merit_value,
                    scaled_product);

    // check convergence or abort if the stepsize becomes too small

    if (quadratic_merit_value
        <= cauchy_quadratic_merit_value + eta * scaled_product)
    {
      (*step_length)                 = alpha;
      (*trial_quadratic_merit_value) = quadratic_merit_value;

      SLEQP_CALL(sleqp_direction_add_scaled(cauchy_direction,
                                            newton_direction,
                                            1. - alpha,
                                            alpha,
                                            zero_eps,
                                            trial_direction));

      break;
    }

    // Update
    {
      alpha *= tau;
    }

    if (alpha <= cutoff_threshold)
    {
      alpha                          = 0.;
      (*step_length)                 = 0.;
      (*trial_quadratic_merit_value) = cauchy_quadratic_merit_value;

      SLEQP_CALL(sleqp_direction_copy(cauchy_direction, trial_direction));

      break;
    }
  }

  SLEQP_CALL(sleqp_timer_stop(linesearch->timer));

#if SLEQP_DEBUG

  {
    bool direction_valid;

    const double zero_eps
      = sleqp_settings_real_value(linesearch->settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

    SLEQP_CALL(sleqp_direction_check(trial_direction,
                                     linesearch->problem,
                                     linesearch->iterate,
                                     multipliers,
                                     linesearch->cache,
                                     zero_eps,
                                     &direction_valid));

    sleqp_num_assert(direction_valid);

    double actual_merit_value;

    SLEQP_CALL(sleqp_merit_quadratic(linesearch->merit,
                                     linesearch->iterate,
                                     trial_direction,
                                     linesearch->penalty_parameter,
                                     &actual_merit_value));

    sleqp_assert_is_eq(actual_merit_value, *trial_quadratic_merit_value, eps);
  }

#endif

  assert(iteration != LINESEARCH_MAX_IT);

  sleqp_log_debug("Cauchy-Newton line search terminated after %d iterations "
                  "(step length: %g, quadratic merit: %g)",
                  iteration,
                  alpha,
                  (*trial_quadratic_merit_value));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
collect_breakpoints(SleqpLineSearch* linesearch,
                    const SleqpVec* violation,
                    const SleqpVec* direction,
                    double* initial_violation,
                    double* initial_slope)
{
  SleqpProblem* problem = linesearch->problem;

  const int num_constraints = sleqp_problem_num_cons(problem);

  const int dim = num_constraints;

  const SleqpVec* v = violation;
  const SleqpVec* d = direction;

  int k_v = 0, k_d = 0;

  while (k_v < v->nnz && k_d < d->nnz)
  {
    bool valid_v = k_v < v->nnz;
    bool valid_d = k_d < d->nnz;

    double d_val = 0., v_val = 0.;

    int idx = valid_v ? v->indices[k_v] : (dim + 1);
    idx     = SLEQP_MIN(idx, valid_d ? d->indices[k_d] : (dim + 1));

    if (valid_v && idx == v->indices[k_v])
    {
      v_val = v->data[k_v++];
    }

    if (valid_d && idx == d->indices[k_d])
    {
      d_val = d->data[k_d++];
    }

    (*initial_violation) += SLEQP_MAX(v_val, 0.);

    if (v_val >= 0.)
    {
      (*initial_slope) += d_val;

      if (d_val < 0.)
      {
        const double breakpoint = (-1.) * v_val / d_val;

        assert(breakpoint >= 0.);

        if (breakpoint <= 1.)
        {
          linesearch->breakpoints[linesearch->num_breakpoints++]
            = (BreakPoint){.point = breakpoint, .slope_change = (-1.) * d_val};
        }
      }
    }
    else
    {
      if (d_val > 0)
      {
        const double breakpoint = (-1.) * v_val / d_val;

        assert(breakpoint >= 0.);

        if (breakpoint <= 1.)
        {
          linesearch->breakpoints[linesearch->num_breakpoints++]
            = (BreakPoint){.point = breakpoint, .slope_change = d_val};
        }
      }
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_breakpoints(SleqpLineSearch* linesearch,
                    double* initial_violation,
                    double* initial_slope)
{
  SleqpProblem* problem = linesearch->problem;

  const double zero_eps
    = sleqp_settings_real_value(linesearch->settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  const int num_constraints = sleqp_problem_num_cons(problem);

  const int max_num_breakpoints = 2 * num_constraints + 1;

  // reserve a sufficient number
  if (linesearch->max_num_breakpoints < max_num_breakpoints)
  {
    SLEQP_CALL(sleqp_realloc(&linesearch->breakpoints, max_num_breakpoints));
    linesearch->max_num_breakpoints = max_num_breakpoints;
  }

  *initial_violation = 0.;
  *initial_slope     = 0.;

  linesearch->num_breakpoints = 0;

  // compute upper violation and break points
  {
    SLEQP_CALL(sleqp_vec_add_scaled(linesearch->cauchy_cons_val,
                                    sleqp_problem_cons_ub(problem),
                                    1.,
                                    -1.,
                                    zero_eps,
                                    linesearch->cauchy_violation));

    SLEQP_CALL(collect_breakpoints(linesearch,
                                   linesearch->cauchy_violation,
                                   linesearch->cauchy_newton_jacobian_prod,
                                   initial_violation,
                                   initial_slope));
  }

  {
    SLEQP_CALL(sleqp_vec_add_scaled(linesearch->cauchy_cons_val,
                                    sleqp_problem_cons_lb(problem),
                                    -1.,
                                    1.,
                                    zero_eps,
                                    linesearch->cauchy_violation));

    SLEQP_CALL(sleqp_vec_scale(linesearch->cauchy_newton_jacobian_prod, -1.));

    SLEQP_CALL(collect_breakpoints(linesearch,
                                   linesearch->cauchy_violation,
                                   linesearch->cauchy_newton_jacobian_prod,
                                   initial_violation,
                                   initial_slope));
  }

  linesearch->breakpoints[linesearch->num_breakpoints++]
    = (BreakPoint){.point = 1., .slope_change = 0.};

  qsort((void*)linesearch->breakpoints,
        linesearch->num_breakpoints,
        sizeof(BreakPoint),
        compare_breakpoints);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_linesearch_trial_step_exact(SleqpLineSearch* linesearch,
                                  const SleqpDirection* cauchy_direction,
                                  const double cauchy_quadratic_merit_value,
                                  const SleqpDirection* newton_direction,
                                  const SleqpVec* multipliers,
                                  SleqpDirection* trial_direction,
                                  double* step_length,
                                  double* trial_quadratic_merit_value)
{
  SleqpIterate* iterate = linesearch->iterate;

  const double zero_eps
    = sleqp_settings_real_value(linesearch->settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  const double eps = sleqp_settings_real_value(linesearch->settings, SLEQP_SETTINGS_REAL_EPS);

  SLEQP_NUM_ASSERT_PARAM(eps);

  const double penalty_parameter = linesearch->penalty_parameter;

  SleqpVec* cauchy_step = sleqp_direction_primal(cauchy_direction);
  SleqpVec* newton_step = sleqp_direction_primal(newton_direction);

  SleqpVec* cauchy_hessian_step = sleqp_direction_hess(cauchy_direction);
  SleqpVec* newton_hessian_step = sleqp_direction_hess(newton_direction);

  SleqpVec* cauchy_cons_jac = sleqp_direction_cons_jac(cauchy_direction);
  SleqpVec* newton_cons_jac = sleqp_direction_cons_jac(newton_direction);

  // Compute Cauchy constraint values
  {
    SLEQP_CALL(sleqp_vec_add(sleqp_iterate_cons_val(iterate),
                             cauchy_cons_jac,
                             zero_eps,
                             linesearch->cauchy_cons_val));
  }

  // Compute Cauchy-Newton (and -Jacobian) direction
  {
    SLEQP_CALL(sleqp_vec_add_scaled(newton_step,
                                    cauchy_step,
                                    1.,
                                    -1.,
                                    zero_eps,
                                    linesearch->cauchy_newton_direction));

    SLEQP_CALL(sleqp_vec_add_scaled(newton_cons_jac,
                                    cauchy_cons_jac,
                                    1.,
                                    -1.,
                                    zero_eps,
                                    linesearch->cauchy_newton_jacobian_prod));
  }

  double offset = 0.;

  // Compute offset
  {
    offset += sleqp_iterate_obj_val(iterate);

    double objective_dot;

    SLEQP_CALL(sleqp_vec_dot(sleqp_iterate_obj_grad(iterate),
                             cauchy_step,
                             &objective_dot));

    offset += objective_dot;

    double cauchy_hessian_dot;

    SLEQP_CALL(
      sleqp_vec_dot(cauchy_step, cauchy_hessian_step, &cauchy_hessian_dot));

    offset += .5 * cauchy_hessian_dot;
  }

  double linear_term = 0.;

  // Compute misc. linear terms
  {
    double objective_dot;

    SLEQP_CALL(sleqp_vec_dot(sleqp_iterate_obj_grad(iterate),
                             linesearch->cauchy_newton_direction,
                             &objective_dot));

    linear_term += objective_dot;

    double mixed_dot;

    SLEQP_CALL(sleqp_vec_dot(linesearch->cauchy_newton_direction,
                             cauchy_hessian_step,
                             &mixed_dot));

    linear_term += mixed_dot;
  }

  double quadratic_term = 0.;

  // Compute quadratic term
  {
    double cauchy_hessian_dot, newton_hessian_dot;

    SLEQP_CALL(sleqp_vec_dot(linesearch->cauchy_newton_direction,
                             cauchy_hessian_step,
                             &cauchy_hessian_dot));

    SLEQP_CALL(sleqp_vec_dot(linesearch->cauchy_newton_direction,
                             newton_hessian_step,
                             &newton_hessian_dot));

    quadratic_term += (newton_hessian_dot - cauchy_hessian_dot);
  }

  // Compute breakpoints
  {
    double violation_value, violation_slope;

    SLEQP_CALL(
      compute_breakpoints(linesearch, &violation_value, &violation_slope));

    offset += penalty_parameter * violation_value;

    linear_term += penalty_parameter * violation_slope;
  }

  double min_value  = offset;
  double min_alpha  = 0.;
  double last_point = 0.;

#if SLEQP_DEBUG
  {
    const double alpha           = last_point;
    const double quadratic_value = offset + .5 * alpha * quadratic_term * alpha;

    SLEQP_NUM_ASSERT_PARAM(quadratic_value);

    double actual_quadratic_merit_value;

    SLEQP_CALL(sleqp_direction_add_scaled(cauchy_direction,
                                          newton_direction,
                                          1. - alpha,
                                          alpha,
                                          zero_eps,
                                          linesearch->test_direction));

    SLEQP_CALL(sleqp_merit_quadratic(linesearch->merit,
                                     iterate,
                                     linesearch->test_direction,
                                     penalty_parameter,
                                     &actual_quadratic_merit_value));

    sleqp_assert_is_eq(actual_quadratic_merit_value, quadratic_value, eps);
  }
#endif

  for (int k = 0; k < linesearch->num_breakpoints; ++k)
  {
    const BreakPoint* breakpoint = linesearch->breakpoints + k;

    const double current_point = breakpoint->point;
    const double point_diff    = current_point - last_point;

    assert(last_point <= current_point);

    // solve quadratic within [last_point, current_point] if we can
    if (quadratic_term != 0.)
    {
      const double alpha = (-1.) * linear_term / quadratic_term;

      const double quadratic_value = offset + linear_term * (alpha - last_point)
                                     + .5 * alpha * quadratic_term * alpha;

      if (last_point <= alpha && alpha <= current_point
          && quadratic_value < min_value)
      {
        min_value = quadratic_value;
        min_alpha = alpha;
      }
    }

    // test current endpoint
    {
      const double alpha = current_point;

      const double quadratic_value = offset + linear_term * (alpha - last_point)
                                     + .5 * alpha * quadratic_term * alpha;

#if SLEQP_DEBUG
      {
        double actual_quadratic_merit_value;

        SLEQP_CALL(sleqp_direction_add_scaled(cauchy_direction,
                                              newton_direction,
                                              1. - alpha,
                                              alpha,
                                              zero_eps,
                                              linesearch->test_direction));

        SLEQP_CALL(sleqp_merit_quadratic(linesearch->merit,
                                         iterate,
                                         linesearch->test_direction,
                                         penalty_parameter,
                                         &actual_quadratic_merit_value));

        sleqp_assert_is_eq(actual_quadratic_merit_value, quadratic_value, eps);
      }
#endif

      if (quadratic_value < min_value)
      {
        min_value = quadratic_value;
        min_alpha = alpha;
      }
    }

    offset += point_diff * linear_term;
    linear_term += penalty_parameter * breakpoint->slope_change;

    last_point = current_point;
  }

  SLEQP_CALL(sleqp_direction_add_scaled(cauchy_direction,
                                        newton_direction,
                                        1. - min_alpha,
                                        min_alpha,
                                        zero_eps,
                                        trial_direction));

  *step_length                 = min_alpha;
  *trial_quadratic_merit_value = min_value;

#if SLEQP_DEBUG

  // Check quadratic merit
  {
    double actual_quadratic_merit;

    SLEQP_CALL(sleqp_merit_quadratic(linesearch->merit,
                                     iterate,
                                     trial_direction,
                                     penalty_parameter,
                                     &actual_quadratic_merit));

    sleqp_assert_is_eq(*(trial_quadratic_merit_value),
                       actual_quadratic_merit,
                       eps);
  }

#endif

  sleqp_log_debug("Cauchy-Newton line search terminated with step length: %g, "
                  "quadratic merit: %g",
                  min_alpha,
                  (*trial_quadratic_merit_value));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
linesearch_free(SleqpLineSearch** star)
{
  SleqpLineSearch* linesearch = *star;

  SLEQP_CALL(sleqp_timer_free(&linesearch->timer));

  sleqp_free(&linesearch->breakpoints);

  SLEQP_CALL(sleqp_direction_release(&linesearch->test_direction));
  sleqp_free(&linesearch->cache);

  SLEQP_CALL(sleqp_vec_free(&linesearch->violated_multipliers));
  SLEQP_CALL(sleqp_vec_free(&linesearch->cauchy_newton_direction));

  SLEQP_CALL(sleqp_vec_free(&linesearch->combined_cons_val));
  SLEQP_CALL(sleqp_vec_free(&linesearch->cauchy_violation));
  SLEQP_CALL(sleqp_vec_free(&linesearch->cauchy_cons_val));
  SLEQP_CALL(sleqp_vec_free(&linesearch->cauchy_newton_jacobian_prod));
  SLEQP_CALL(sleqp_vec_free(&linesearch->cauchy_point));

  SLEQP_CALL(sleqp_merit_release(&linesearch->merit));

  SLEQP_CALL(sleqp_iterate_release(&linesearch->iterate));

  SLEQP_CALL(sleqp_settings_release(&linesearch->settings));

  SLEQP_CALL(sleqp_problem_release(&linesearch->problem));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SleqpTimer*
sleqp_linesearch_get_timer(SleqpLineSearch* linesearch)
{
  return linesearch->timer;
}

SLEQP_RETCODE
sleqp_linesearch_capture(SleqpLineSearch* linesearch)
{
  ++linesearch->refcount;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_linesearch_release(SleqpLineSearch** star)
{
  SleqpLineSearch* linesearch = *star;

  if (!linesearch)
  {
    return SLEQP_OKAY;
  }

  if (--linesearch->refcount == 0)
  {
    SLEQP_CALL(linesearch_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
