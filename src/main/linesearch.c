#include "linesearch.h"

#include "cmp.h"
#include "fail.h"
#include "feas.h"
#include "mem.h"
#include "util.h"

#define LINESEARCH_MAX_IT 10000

typedef struct BreakPoint
{
  double slope_change;
  double point;
} BreakPoint;

static int compare_breakpoints(const void* first, const void* second)
{
  const double first_point = ((const BreakPoint*) first)->point;
  const double second_point = ((const BreakPoint*) second)->point;

  if(first_point < second_point)
  {
    return -1;
  }
  else if(first_point > second_point)
  {
    return 1;
  }

  return 0;
}

struct SleqpLineSearchData
{
  int refcount;

  SleqpProblem* problem;
  SleqpParams* params;
  SleqpMerit* merit;

  SleqpIterate* iterate;
  double penalty_parameter;
  double trust_radius;

  double* prod_cache;

  SleqpSparseVec* cauchy_point;

  SleqpSparseVec* cauchy_jacobian_prod;
  SleqpSparseVec* cauchy_newton_jacobian_prod;
  SleqpSparseVec* newton_jacobian_prod;
  SleqpSparseVec* cauchy_cons_val;

  SleqpSparseVec* cauchy_violation;

  SleqpSparseVec* combined_cons_val;

  SleqpSparseVec* cauchy_newton_direction;

  SleqpSparseVec* violated_multipliers;

  SleqpSparseVec* test_direction;

  int num_breakpoints;
  int max_num_breakpoints;
  BreakPoint* breakpoints;

  SleqpTimer* timer;
};

SLEQP_RETCODE sleqp_linesearch_create(SleqpLineSearchData** star,
                                      SleqpProblem* problem,
                                      SleqpParams* params,
                                      SleqpMerit* merit)
{
  SLEQP_CALL(sleqp_malloc(star));

  const int num_variables = sleqp_problem_num_variables(problem);

  const int num_constraints = sleqp_problem_num_constraints(problem);

  SleqpLineSearchData* linesearch = *star;

  *linesearch = (SleqpLineSearchData){0};

  linesearch->refcount = 1;

  linesearch->problem = problem;
  SLEQP_CALL(sleqp_problem_capture(linesearch->problem));

  SLEQP_CALL(sleqp_params_capture(params));
  linesearch->params = params;

  SLEQP_CALL(sleqp_merit_capture(merit));

  linesearch->merit = merit;

  SLEQP_CALL(sleqp_alloc_array(&linesearch->prod_cache,
                               num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&linesearch->cauchy_point,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&linesearch->cauchy_jacobian_prod,
                                              num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&linesearch->cauchy_newton_jacobian_prod,
                                              num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&linesearch->newton_jacobian_prod,
                                              num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&linesearch->cauchy_cons_val,
                                              num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&linesearch->cauchy_violation,
                                              num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&linesearch->combined_cons_val,
                                              num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&linesearch->cauchy_newton_direction,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&linesearch->violated_multipliers,
                                              num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&linesearch->test_direction,
                                              num_variables));

  SLEQP_CALL(sleqp_timer_create(&linesearch->timer));

  linesearch->max_num_breakpoints = 0;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_linesearch_set_iterate(SleqpLineSearchData* linesearch,
                                           SleqpIterate* iterate,
                                           double penalty_parameter,
                                           double trust_radius)
{
  SLEQP_CALL(sleqp_iterate_release(&linesearch->iterate));

  linesearch->iterate = iterate;

  SLEQP_CALL(sleqp_iterate_capture(linesearch->iterate));

  linesearch->penalty_parameter = penalty_parameter;
  linesearch->trust_radius = trust_radius;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_linesearch_cauchy_step(SleqpLineSearchData* linesearch,
                                           SleqpSparseVec* direction,
                                           const SleqpSparseVec* multipliers,
                                           SleqpSparseVec* hessian_direction,
                                           bool* full_step,
                                           double* quadratic_merit_value)
{
  SleqpProblem* problem = linesearch->problem;
  SleqpMerit* merit_data = linesearch->merit;

  const int num_constraints = sleqp_problem_num_constraints(problem);

  SleqpIterate* iterate = linesearch->iterate;
  const double penalty_parameter = linesearch->penalty_parameter;
  const double trust_radius = linesearch->trust_radius;

  const double one = 1.;

  const double eps = sleqp_params_get(linesearch->params,
                                      SLEQP_PARAM_EPS);

  const double zero_eps = sleqp_params_get(linesearch->params,
                                           SLEQP_PARAM_ZERO_EPS);

  SLEQP_CALL(sleqp_timer_start(linesearch->timer));

  double exact_violation;

  SLEQP_CALL(sleqp_violation_one_norm(problem,
                                      sleqp_iterate_get_cons_val(iterate),
                                      &exact_violation));

#if !defined(NDEBUG)

  // Check Hessian product
  {
    SLEQP_CALL(sleqp_problem_hess_prod(problem,
                                       &one,
                                       direction,
                                       multipliers,
                                       linesearch->test_direction));

    sleqp_num_assert(sleqp_sparse_vector_eq(hessian_direction,
                                            linesearch->test_direction,
                                            eps));

  }

#endif

  (*full_step) = true;
  (*quadratic_merit_value) = 0.;

  const double exact_merit_value = sleqp_iterate_get_func_val(iterate) + penalty_parameter*exact_violation;

  double hessian_product;

  SLEQP_CALL(sleqp_sparse_vector_dot(direction,
                                     hessian_direction,
                                     &hessian_product));

  double objective_dot;
  SleqpSparseVec* jacobian_product = linesearch->cauchy_jacobian_prod;

  double delta = 1.;

  {
    double direction_norm = sleqp_sparse_vector_norm(direction);

    double direction_factor = trust_radius / direction_norm;

    if(direction_factor < 1)
    {
      delta = direction_factor;

      hessian_product *= (delta*delta);

      SLEQP_CALL(sleqp_sparse_vector_scale(direction, delta));
    }
  }

  // prepare initial products
  {
    SLEQP_CALL(sleqp_sparse_vector_dot(sleqp_iterate_get_func_grad(iterate),
                                       direction,
                                       &objective_dot));

    SLEQP_CALL(sleqp_sparse_matrix_vector_product(sleqp_iterate_get_cons_jac(iterate),
                                                  direction,
                                                  linesearch->prod_cache));

    SLEQP_CALL(sleqp_sparse_vector_from_raw(jacobian_product,
                                            linesearch->prod_cache,
                                            num_constraints,
                                            zero_eps));
  }

  const double eta = sleqp_params_get(linesearch->params, SLEQP_PARAM_CAUCHY_ETA);
  const double tau = sleqp_params_get(linesearch->params, SLEQP_PARAM_CAUCHY_TAU);

  double linear_violation;

  int iteration = 0;

  sleqp_log_debug("Beginning Cauchy line search, exact merit value: %f",
                  exact_merit_value);

  for(iteration = 0; iteration < LINESEARCH_MAX_IT; ++iteration)
  {
    double linear_merit_value = 0.;

    // compute linear merit value
    {
      linear_merit_value += sleqp_iterate_get_func_val(iterate) + objective_dot;

      SLEQP_CALL(sleqp_sparse_vector_add(sleqp_iterate_get_cons_val(iterate),
                                         jacobian_product,
                                         zero_eps,
                                         linesearch->combined_cons_val));

      SLEQP_CALL(sleqp_violation_one_norm(problem,
                                          linesearch->combined_cons_val,
                                          &linear_violation));

      linear_merit_value += penalty_parameter * linear_violation;
    }

    // compute quadratic merit value
    {
      (*quadratic_merit_value) = linear_merit_value + (0.5 * hessian_product);
    }

    sleqp_log_debug("Cauchy line search iteration %d, step length: %g, linear merit value: %f, quadratic merit value: %f",
                    iteration,
                    delta,
                    linear_merit_value,
                    *quadratic_merit_value);

#if !defined(NDEBUG)

    {
      {
        double actual_linear_merit_value;

        SLEQP_CALL(sleqp_merit_linear(merit_data,
                                      iterate,
                                      direction,
                                      penalty_parameter,
                                      &actual_linear_merit_value));

        sleqp_assert_is_eq(linear_merit_value,
                           actual_linear_merit_value,
                           eps);
      }

      {
        double actual_quadratic_merit_value;

        SLEQP_CALL(sleqp_merit_quadratic(merit_data,
                                         iterate,
                                         &one,
                                         direction,
                                         multipliers,
                                         penalty_parameter,
                                         &actual_quadratic_merit_value));

        sleqp_assert_is_eq((*quadratic_merit_value),
                           actual_quadratic_merit_value,
                           eps);
      }
    }

#endif

    // check condition
    // This is a numerically more sensible version of the condition
    // ((exact_merit_value - (*quadratic_merit_value)) >= eta*(exact_merit_value - linear_merit_value))
    if((penalty_parameter*(exact_violation - linear_violation) - objective_dot) * (1. - eta) >= (0.5 * hessian_product))
    {
      break;
    }

    // update products
    SLEQP_CALL(sleqp_sparse_vector_scale(direction, tau));

    hessian_product *= (tau*tau);

    objective_dot *= tau;

    SLEQP_CALL(sleqp_sparse_vector_scale(jacobian_product, tau));

    delta *= tau;

    if(sleqp_is_zero(delta, eps))
    {
      SLEQP_CALL(sleqp_sparse_vector_clear(direction));
      delta = 0.;
      break;
    }
  }

  SLEQP_CALL(sleqp_sparse_vector_scale(hessian_direction, delta));

  assert(iteration != LINESEARCH_MAX_IT);

  sleqp_log_debug("Cauchy line search terminated after %d iterations (step length: %f, quadratic merit: %f)",
                  iteration,
                  delta,
                  (*quadratic_merit_value));

  (*full_step) = false;

  sleqp_assert_is_leq(sleqp_sparse_vector_norm(direction), trust_radius, eps);

  SLEQP_CALL(sleqp_timer_stop(linesearch->timer));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_linesearch_trial_step(SleqpLineSearchData* linesearch,
                                          const SleqpSparseVec* cauchy_step,
                                          const SleqpSparseVec* cauchy_hessian_step,
                                          const double cauchy_quadratic_merit_value,
                                          const SleqpSparseVec* newton_step,
                                          const SleqpSparseVec* newton_hessian_step,
                                          const SleqpSparseVec* multipliers,
                                          SleqpSparseVec* trial_step,
                                          double* step_length,
                                          double* trial_quadratic_merit_value)
{
  SleqpProblem* problem = linesearch->problem;
  SleqpMerit* merit_data = linesearch->merit;
  SleqpIterate* iterate = linesearch->iterate;
  const double penalty_parameter = linesearch->penalty_parameter;

  const int num_constraints = sleqp_problem_num_constraints(problem);

  const double eps = sleqp_params_get(linesearch->params,
                                      SLEQP_PARAM_EPS);

  SLEQP_NUM_ASSERT_PARAM(eps);

  const double zero_eps = sleqp_params_get(linesearch->params,
                                           SLEQP_PARAM_ZERO_EPS);

  SLEQP_CALL(sleqp_timer_start(linesearch->timer));

  // Compute Cauchy-Newton direction
  {
    SLEQP_CALL(sleqp_sparse_vector_add_scaled(newton_step,
                                              cauchy_step,
                                              1.,
                                              -1.,
                                              zero_eps,
                                              linesearch->cauchy_newton_direction));
  }

#if !defined(NDEBUG)

  // Check merit and hessian products
  {
    double func_dual = 1.;

    double actual_quadratic_merit;

    SLEQP_CALL(sleqp_merit_quadratic(merit_data,
                                     iterate,
                                     &func_dual,
                                     cauchy_step,
                                     multipliers,
                                     penalty_parameter,
                                     &actual_quadratic_merit));

    sleqp_assert_is_eq(cauchy_quadratic_merit_value,
                       actual_quadratic_merit,
                       eps);

    SLEQP_CALL(sleqp_problem_hess_prod(problem,
                                       &func_dual,
                                       cauchy_step,
                                       multipliers,
                                       linesearch->test_direction));

    sleqp_num_assert(sleqp_sparse_vector_eq(cauchy_hessian_step,
                                            linesearch->test_direction,
                                            eps));

    SLEQP_CALL(sleqp_problem_hess_prod(problem,
                                       &func_dual,
                                       newton_step,
                                       multipliers,
                                       linesearch->test_direction));

    sleqp_num_assert(sleqp_sparse_vector_eq(newton_hessian_step,
                                            linesearch->test_direction,
                                            eps));
  }

#endif

  // Compute Jacobian products
  {
    SLEQP_CALL(sleqp_sparse_matrix_vector_product(sleqp_iterate_get_cons_jac(iterate),
                                                  cauchy_step,
                                                  linesearch->prod_cache));

    SLEQP_CALL(sleqp_sparse_vector_from_raw(linesearch->cauchy_jacobian_prod,
                                            linesearch->prod_cache,
                                            num_constraints,
                                            zero_eps));

    SLEQP_CALL(sleqp_sparse_matrix_vector_product(sleqp_iterate_get_cons_jac(iterate),
                                                  newton_step,
                                                  linesearch->prod_cache));

    SLEQP_CALL(sleqp_sparse_vector_from_raw(linesearch->newton_jacobian_prod,
                                            linesearch->prod_cache,
                                            num_constraints,
                                            zero_eps));
  }

  // Compute initial scalar products for the linear merit function
  double cauchy_gradient_dot;

  SLEQP_CALL(sleqp_sparse_vector_dot(sleqp_iterate_get_func_grad(iterate),
                                     cauchy_step,
                                     &cauchy_gradient_dot));

  double newton_gradient_dot;

  SLEQP_CALL(sleqp_sparse_vector_dot(sleqp_iterate_get_func_grad(iterate),
                                     newton_step,
                                     &newton_gradient_dot));

  // Compute initial scalar products for the quadratic merit function
  double cauchy_cauchy_product;

  SLEQP_CALL(sleqp_sparse_vector_dot(cauchy_step,
                                     cauchy_hessian_step,
                                     &cauchy_cauchy_product));

  double cauchy_newton_product;

  SLEQP_CALL(sleqp_sparse_vector_dot(cauchy_step,
                                     newton_hessian_step,
                                     &cauchy_newton_product));

  double newton_newton_product;

  SLEQP_CALL(sleqp_sparse_vector_dot(newton_step,
                                     newton_hessian_step,
                                     &newton_newton_product));

  const double eta = sleqp_params_get(linesearch->params,
                                      SLEQP_PARAM_LINESEARCH_ETA);

  const double tau = sleqp_params_get(linesearch->params,
                                      SLEQP_PARAM_LINESEARCH_TAU);

  const double cutoff_threshold = sleqp_params_get(linesearch->params,
                                                   SLEQP_PARAM_LINESEARCH_CUTOFF);

  int iteration = 0;

  double alpha = 1.;

  // compute the maximum trial step length
  {
    SLEQP_CALL(sleqp_sparse_vector_add(sleqp_iterate_get_primal(iterate),
                                       cauchy_step,
                                       zero_eps,
                                       linesearch->cauchy_point));

    SLEQP_CALL(sleqp_max_step_length(linesearch->cauchy_point,
                                     linesearch->cauchy_newton_direction,
                                     sleqp_problem_var_lb(problem),
                                     sleqp_problem_var_ub(problem),
                                     &alpha));

    assert(alpha >= 0);
    assert(alpha <= 1);
  }

  double quadratic_merit_gradient_product = 0.;

  // compute merit gradient products
  {
    SLEQP_CALL(sleqp_sparse_vector_add(sleqp_iterate_get_cons_val(iterate),
                                       linesearch->cauchy_jacobian_prod,
                                       zero_eps,
                                       linesearch->cauchy_cons_val));

    // use the violated multipliers in 0, +/-1 to filter the Jacobian product
    SLEQP_CALL(sleqp_violated_constraint_multipliers(problem,
                                                     linesearch->cauchy_cons_val,
                                                     linesearch->violated_multipliers,
                                                     NULL));

    double jacobian_dot;

    SLEQP_CALL(sleqp_sparse_vector_dot(linesearch->violated_multipliers,
                                       linesearch->cauchy_jacobian_prod,
                                       &jacobian_dot));

    double quadratic_merit_gradient_cauchy = cauchy_gradient_dot + jacobian_dot + cauchy_cauchy_product;

    SLEQP_CALL(sleqp_sparse_vector_dot(linesearch->violated_multipliers,
                                       linesearch->newton_jacobian_prod,
                                       &jacobian_dot));

    double quadratic_merit_gradient_newton = newton_gradient_dot + jacobian_dot + cauchy_newton_product;

    quadratic_merit_gradient_product = (quadratic_merit_gradient_newton - quadratic_merit_gradient_cauchy);
  }

  if(alpha <= cutoff_threshold)
  {
    (*step_length) = 0.;
    (*trial_quadratic_merit_value) = cauchy_quadratic_merit_value;

    SLEQP_CALL(sleqp_sparse_vector_copy(cauchy_step, trial_step));

    SLEQP_CALL(sleqp_timer_stop(linesearch->timer));

    return SLEQP_OKAY;
  }

  sleqp_log_debug("Beginning Cauchy-Newton line search, initial step length: %g, quadratic merit at Cauchy point: %g",
                  alpha,
                  cauchy_quadratic_merit_value);

  for(iteration = 0; iteration < LINESEARCH_MAX_IT; ++iteration)
  {
    double linear_merit_value = 0.;
    double quadratic_merit_value = 0.;

    // Compute linear merit
    {
      linear_merit_value += sleqp_iterate_get_func_val(iterate)
        + ((1. - alpha) * cauchy_gradient_dot)
        + ((alpha) * newton_gradient_dot);

      SLEQP_CALL(sleqp_sparse_vector_add_scaled(sleqp_iterate_get_cons_val(iterate),
                                                linesearch->cauchy_jacobian_prod,
                                                1.,
                                                (1. - alpha),
                                                zero_eps,
                                                linesearch->cauchy_cons_val));

      SLEQP_CALL(sleqp_sparse_vector_add_scaled(linesearch->cauchy_cons_val,
                                                linesearch->newton_jacobian_prod,
                                                1.,
                                                alpha,
                                                zero_eps,
                                                linesearch->combined_cons_val));

      double total_violation;

      SLEQP_CALL(sleqp_violation_one_norm(problem,
                                          linesearch->combined_cons_val,
                                          &total_violation));

      linear_merit_value += penalty_parameter * total_violation;
    }

    // Compute quadratic term

    {
      double quadratic_term = 0.;

      quadratic_term += 0.5 *(1. - alpha) * (1. - alpha) * cauchy_cauchy_product;
      quadratic_term += alpha * ((1. - alpha) * cauchy_newton_product + 0.5 * alpha * newton_newton_product);

      quadratic_merit_value = linear_merit_value + quadratic_term;
    }

#if !defined(NDEBUG)

    {
      SLEQP_CALL(sleqp_sparse_vector_add_scaled(cauchy_step,
                                                newton_step,
                                                1. - alpha,
                                                alpha,
                                                zero_eps,
                                                linesearch->test_direction));
    }

    {
      double actual_linear_merit;

      SLEQP_CALL(sleqp_merit_linear(merit_data,
                                    iterate,
                                    linesearch->test_direction,
                                    penalty_parameter,
                                    &actual_linear_merit));

      sleqp_assert_is_eq(linear_merit_value,
                         actual_linear_merit,
                         eps);
    }

    // Check quadratic merit
    {
      double one = 1.;

      double actual_quadratic_merit;

      SLEQP_CALL(sleqp_merit_quadratic(merit_data,
                                       iterate,
                                       &one,
                                       linesearch->test_direction,
                                       multipliers,
                                       penalty_parameter,
                                       &actual_quadratic_merit));

      sleqp_assert_is_eq(quadratic_merit_value,
                         actual_quadratic_merit,
                         eps);
    }

#endif

    const double scaled_product = alpha*quadratic_merit_gradient_product;

    sleqp_log_debug("Cauchy-Newton line search iteration %d, step length: %g, quadratic merit value: %g, scaled inner product: %g",
                    iteration,
                    alpha,
                    quadratic_merit_value,
                    scaled_product);

    // check convergence or abort if the stepsize becomes too small

    if(quadratic_merit_value <= cauchy_quadratic_merit_value + eta*scaled_product)
    {
      (*step_length) = alpha;
      (*trial_quadratic_merit_value) = quadratic_merit_value;

      SLEQP_CALL(sleqp_sparse_vector_add_scaled(cauchy_step,
                                                newton_step,
                                                1. - alpha,
                                                alpha,
                                                zero_eps,
                                                trial_step));

      break;
    }

    // Update
    {
      alpha *= tau;
    }

    if(alpha <= cutoff_threshold)
    {
      alpha = 0.;
      (*step_length) = 0.;
      (*trial_quadratic_merit_value) = cauchy_quadratic_merit_value;

      SLEQP_CALL(sleqp_sparse_vector_copy(cauchy_step, trial_step));

      break;
    }

  }

  SLEQP_CALL(sleqp_timer_stop(linesearch->timer));

  assert(iteration != LINESEARCH_MAX_IT);

  sleqp_log_debug("Cauchy-Newton line search terminated after %d iterations (step length: %g, quadratic merit: %g)",
                  iteration,
                  alpha,
                  (*trial_quadratic_merit_value));

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE collect_breakpoints(SleqpLineSearchData* linesearch,
                                  const SleqpSparseVec* violation,
                                  const SleqpSparseVec* direction,
                                  double* initial_violation,
                                  double* initial_slope)
{
  SleqpProblem* problem = linesearch->problem;

  const int num_constraints = sleqp_problem_num_constraints(problem);

  const int dim = num_constraints;

  const SleqpSparseVec* v = violation;
  const SleqpSparseVec* d = direction;

  int k_v = 0, k_d = 0;

  while(k_v < v->nnz && k_d < d->nnz)
  {
    bool valid_v = k_v < v->nnz;
    bool valid_d = k_d < d->nnz;

    double d_val = 0., v_val = 0.;

    int idx = valid_v ? v->indices[k_v] : (dim + 1);
    idx = SLEQP_MIN(idx, valid_d ? d->indices[k_d] : (dim + 1));

    if(valid_v && idx == v->indices[k_v])
    {
      v_val = v->data[k_v++];
    }

    if(valid_d && idx == d->indices[k_d])
    {
      d_val = d->data[k_d++];
    }

    (*initial_violation) += SLEQP_MAX(v_val, 0.);

    if(v_val >= 0.)
    {
      (*initial_slope) += d_val;

      if(d_val < 0.)
      {
        const double breakpoint = (-1.) * v_val / d_val;

        assert(breakpoint >= 0.);

        if(breakpoint <= 1.)
        {
          linesearch->breakpoints[linesearch->num_breakpoints++] =
            (BreakPoint) { .point = breakpoint, .slope_change = (-1.) * d_val };
        }
      }
    }
    else
    {
      if(d_val > 0)
      {
        const double breakpoint = (-1.) * v_val / d_val;

        assert(breakpoint >= 0.);

        if(breakpoint <= 1.)
        {
          linesearch->breakpoints[linesearch->num_breakpoints++] =
            (BreakPoint) { .point = breakpoint, .slope_change = d_val };
        }
      }
    }
  }

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE compute_breakpoints(SleqpLineSearchData* linesearch,
                                  double* initial_violation,
                                  double* initial_slope)
{
  SleqpProblem* problem = linesearch->problem;

  const double zero_eps = sleqp_params_get(linesearch->params,
                                           SLEQP_PARAM_ZERO_EPS);

  const int num_constraints = sleqp_problem_num_constraints(problem);

  const int max_num_breakpoints = 2*num_constraints + 1;

  // reserve a sufficient number
  if(linesearch->max_num_breakpoints < max_num_breakpoints)
  {
    SLEQP_CALL(sleqp_realloc(&linesearch->breakpoints, max_num_breakpoints));
    linesearch->max_num_breakpoints = max_num_breakpoints;
  }

  *initial_violation = 0.;
  *initial_slope = 0.;

  linesearch->num_breakpoints = 0;

  // compute upper violation and break points
  {
    SLEQP_CALL(sleqp_sparse_vector_add_scaled(linesearch->cauchy_cons_val,
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
    SLEQP_CALL(sleqp_sparse_vector_add_scaled(linesearch->cauchy_cons_val,
                                              sleqp_problem_cons_lb(problem),
                                              -1.,
                                              1.,
                                              zero_eps,
                                              linesearch->cauchy_violation));

    SLEQP_CALL(sleqp_sparse_vector_scale(linesearch->cauchy_newton_jacobian_prod,
                                         -1.));

    SLEQP_CALL(collect_breakpoints(linesearch,
                                   linesearch->cauchy_violation,
                                   linesearch->cauchy_newton_jacobian_prod,
                                   initial_violation,
                                   initial_slope));
  }

  linesearch->breakpoints[linesearch->num_breakpoints++] =
    (BreakPoint) { .point = 1., .slope_change = 0. };

  qsort((void*) linesearch->breakpoints,
        linesearch->num_breakpoints,
        sizeof(BreakPoint),
        compare_breakpoints);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_linesearch_trial_step_exact(SleqpLineSearchData* linesearch,
                                                const SleqpSparseVec* cauchy_step,
                                                const SleqpSparseVec* cauchy_hessian_step,
                                                const double cauchy_quadratic_merit_value,
                                                const SleqpSparseVec* newton_step,
                                                const SleqpSparseVec* newton_hessian_step,
                                                const SleqpSparseVec* multipliers,
                                                SleqpSparseVec* trial_step,
                                                double* step_length,
                                                double* trial_quadratic_merit_value)
{
  SleqpProblem* problem = linesearch->problem;
  SleqpIterate* iterate = linesearch->iterate;

  const double zero_eps = sleqp_params_get(linesearch->params,
                                           SLEQP_PARAM_ZERO_EPS);

  const double eps = sleqp_params_get(linesearch->params,
                                      SLEQP_PARAM_EPS);

  SLEQP_NUM_ASSERT_PARAM(eps);

  const double penalty_parameter = linesearch->penalty_parameter;

  const int num_constraints = sleqp_problem_num_constraints(problem);

  // Compute Cauchy constraint values
  {
    SLEQP_CALL(sleqp_sparse_matrix_vector_product(sleqp_iterate_get_cons_jac(iterate),
                                                  cauchy_step,
                                                  linesearch->prod_cache));

    SLEQP_CALL(sleqp_sparse_vector_from_raw(linesearch->cauchy_jacobian_prod,
                                            linesearch->prod_cache,
                                            num_constraints,
                                            zero_eps));

    SLEQP_CALL(sleqp_sparse_vector_add(sleqp_iterate_get_cons_val(iterate),
                                       linesearch->cauchy_jacobian_prod,
                                       zero_eps,
                                       linesearch->cauchy_cons_val));
  }

  // Compute Cauchy-Newton (and -Jacobian) direction
  {
    SLEQP_CALL(sleqp_sparse_vector_add_scaled(newton_step,
                                              cauchy_step,
                                              1.,
                                              -1.,
                                              zero_eps,
                                              linesearch->cauchy_newton_direction));

    SLEQP_CALL(sleqp_sparse_matrix_vector_product(sleqp_iterate_get_cons_jac(iterate),
                                                  linesearch->cauchy_newton_direction,
                                                  linesearch->prod_cache));

    SLEQP_CALL(sleqp_sparse_vector_from_raw(linesearch->cauchy_newton_jacobian_prod,
                                            linesearch->prod_cache,
                                            num_constraints,
                                            zero_eps));
  }

  double offset = 0.;

  // Compute offset
  {
    offset += sleqp_iterate_get_func_val(iterate);

    double objective_dot;

    SLEQP_CALL(sleqp_sparse_vector_dot(sleqp_iterate_get_func_grad(iterate),
                                       cauchy_step,
                                       &objective_dot));

    offset += objective_dot;

    double cauchy_hessian_dot;

    SLEQP_CALL(sleqp_sparse_vector_dot(cauchy_step, cauchy_hessian_step, &cauchy_hessian_dot));

    offset += .5 * cauchy_hessian_dot;
  }

  double linear_term = 0.;

  // Compute misc. linear terms
  {
    double objective_dot;

    SLEQP_CALL(sleqp_sparse_vector_dot(sleqp_iterate_get_func_grad(iterate),
                                       linesearch->cauchy_newton_direction,
                                       &objective_dot));

    linear_term += objective_dot;

    double mixed_dot;

    SLEQP_CALL(sleqp_sparse_vector_dot(linesearch->cauchy_newton_direction,
                                       cauchy_hessian_step,
                                       &mixed_dot));

    linear_term += mixed_dot;
  }

  double quadratic_term = 0.;

  // Compute quadratic term
  {
    double cauchy_hessian_dot, newton_hessian_dot;

    SLEQP_CALL(sleqp_sparse_vector_dot(linesearch->cauchy_newton_direction,
                                       cauchy_hessian_step,
                                       &cauchy_hessian_dot));

    SLEQP_CALL(sleqp_sparse_vector_dot(linesearch->cauchy_newton_direction,
                                       newton_hessian_step,
                                       &newton_hessian_dot));

    quadratic_term += (newton_hessian_dot - cauchy_hessian_dot);
  }

  // Compute breakpoints
  {
    double violation_value, violation_slope;

    SLEQP_CALL(compute_breakpoints(linesearch,
                                   &violation_value,
                                   &violation_slope));

    offset += penalty_parameter * violation_value;

    linear_term += penalty_parameter * violation_slope;
  }

  double min_value = offset;
  double min_alpha = 0.;
  double last_point = 0.;

#if !defined(NDEBUG)
  {
    const double one = 1.;
    const double alpha = last_point;
    const double quadratic_value = offset + .5 * alpha * quadratic_term * alpha;

    SLEQP_NUM_ASSERT_PARAM(quadratic_value);

    double actual_quadratic_merit_value;

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(cauchy_step,
                                              newton_step,
                                              1. - alpha,
                                              alpha,
                                              zero_eps,
                                              linesearch->test_direction));

    SLEQP_CALL(sleqp_merit_quadratic(linesearch->merit,
                                     iterate,
                                     &one,
                                     linesearch->test_direction,
                                     multipliers,
                                     penalty_parameter,
                                     &actual_quadratic_merit_value));

    sleqp_assert_is_eq(actual_quadratic_merit_value, quadratic_value, eps);
  }
#endif

  for(int k = 0; k < linesearch->num_breakpoints; ++k)
  {
    const BreakPoint* breakpoint = linesearch->breakpoints + k;

    const double current_point = breakpoint->point;
    const double point_diff = current_point - last_point;

    assert(last_point <= current_point);

    // solve quadratic within [last_point, current_point] if we can
    if(quadratic_term != 0.)
    {
      const double alpha = (-1.) * linear_term / quadratic_term;

      const double quadratic_value = offset
        + linear_term * (alpha - last_point)
        + .5 * alpha * quadratic_term * alpha;

      if(last_point <= alpha &&
         alpha <= current_point &&
         quadratic_value < min_value)
      {
        min_value = quadratic_value;
        min_alpha = alpha;
      }
    }

    // test current endpoint
    {
      const double alpha = current_point;

      const double quadratic_value = offset
        + linear_term * (alpha - last_point)
        + .5 * alpha * quadratic_term * alpha;

#if !defined(NDEBUG)
      {
        const double one = 1.;

        double actual_quadratic_merit_value;

        SLEQP_CALL(sleqp_sparse_vector_add_scaled(cauchy_step,
                                                  newton_step,
                                                  1. - alpha,
                                                  alpha,
                                                  zero_eps,
                                                  linesearch->test_direction));

        SLEQP_CALL(sleqp_merit_quadratic(linesearch->merit,
                                         iterate,
                                         &one,
                                         linesearch->test_direction,
                                         multipliers,
                                         penalty_parameter,
                                         &actual_quadratic_merit_value));

        sleqp_assert_is_eq(actual_quadratic_merit_value, quadratic_value, eps);
      }
#endif

      if(quadratic_value < min_value)
      {
        min_value = quadratic_value;
        min_alpha = alpha;
      }
    }

    offset += point_diff * linear_term;
    linear_term += penalty_parameter * breakpoint->slope_change;

    last_point = current_point;
  }

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(cauchy_step,
                                            newton_step,
                                            1. - min_alpha,
                                            min_alpha,
                                            zero_eps,
                                            trial_step));

  *step_length = min_alpha;
  *trial_quadratic_merit_value = min_value;

#if !defined(NDEBUG)

  // Check quadratic merit
  {
    double one = 1.;

    double actual_quadratic_merit;

    SLEQP_CALL(sleqp_merit_quadratic(linesearch->merit,
                                     iterate,
                                     &one,
                                     trial_step,
                                     multipliers,
                                     penalty_parameter,
                                     &actual_quadratic_merit));

    sleqp_assert_is_eq(*(trial_quadratic_merit_value),
                       actual_quadratic_merit,
                       eps);
  }

#endif

  sleqp_log_debug("Cauchy-Newton line search terminated with step length: %g, quadratic merit: %g",
                  min_alpha,
                  (*trial_quadratic_merit_value));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE linesearch_free(SleqpLineSearchData** star)
{
  SleqpLineSearchData* linesearch = *star;

  SLEQP_CALL(sleqp_timer_free(&linesearch->timer));

  sleqp_free(&linesearch->breakpoints);

  SLEQP_CALL(sleqp_sparse_vector_free(&linesearch->test_direction));

  SLEQP_CALL(sleqp_sparse_vector_free(&linesearch->violated_multipliers));
  SLEQP_CALL(sleqp_sparse_vector_free(&linesearch->cauchy_newton_direction));

  SLEQP_CALL(sleqp_sparse_vector_free(&linesearch->combined_cons_val));
  SLEQP_CALL(sleqp_sparse_vector_free(&linesearch->cauchy_violation));
  SLEQP_CALL(sleqp_sparse_vector_free(&linesearch->cauchy_cons_val));
  SLEQP_CALL(sleqp_sparse_vector_free(&linesearch->newton_jacobian_prod));
  SLEQP_CALL(sleqp_sparse_vector_free(&linesearch->cauchy_newton_jacobian_prod));
  SLEQP_CALL(sleqp_sparse_vector_free(&linesearch->cauchy_jacobian_prod));
  SLEQP_CALL(sleqp_sparse_vector_free(&linesearch->cauchy_point));
  sleqp_free(&linesearch->prod_cache);

  SLEQP_CALL(sleqp_merit_release(&linesearch->merit));

  SLEQP_CALL(sleqp_iterate_release(&linesearch->iterate));

  SLEQP_CALL(sleqp_params_release(&linesearch->params));

  SLEQP_CALL(sleqp_problem_release(&linesearch->problem));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SleqpTimer* sleqp_linesearch_get_timer(SleqpLineSearchData* linesearch)
{
  return linesearch->timer;
}

SLEQP_RETCODE sleqp_linesearch_capture(SleqpLineSearchData* linesearch)
{
  ++linesearch->refcount;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_linesearch_release(SleqpLineSearchData** star)
{
  SleqpLineSearchData* linesearch = *star;

  if(!linesearch)
  {
    return SLEQP_OKAY;
  }

  if(--linesearch->refcount == 0)
  {
    SLEQP_CALL(linesearch_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
