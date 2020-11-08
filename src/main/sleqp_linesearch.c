#include "sleqp_linesearch.h"

#include "sleqp_assert.h"
#include "sleqp_feas.h"

#define LINESEARCH_MAX_IT 10000

struct SleqpLineSearchData
{
  int refcount;

  SleqpProblem* problem;
  SleqpParams* params;
  SleqpMeritData* merit_data;

  SleqpIterate* iterate;
  double penalty_parameter;
  double trust_radius;

  double* prod_cache;

  SleqpSparseVec* cauchy_point;

  SleqpSparseVec* cauchy_jacobian_prod;
  SleqpSparseVec* newton_jacobian_prod;
  SleqpSparseVec* cauchy_cons_val;

  SleqpSparseVec* combined_cons_val;

  SleqpSparseVec* cauchy_newton_direction;

  SleqpSparseVec* violated_multipliers;

  SleqpSparseVec* test_direction;
};

SLEQP_RETCODE sleqp_linesearch_create(SleqpLineSearchData** star,
                                      SleqpProblem* problem,
                                      SleqpParams* params,
                                      SleqpMeritData* merit_data)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpLineSearchData* linesearch = *star;

  *linesearch = (SleqpLineSearchData){0};

  linesearch->refcount = 1;

  linesearch->problem = problem;

  SLEQP_CALL(sleqp_params_capture(params));
  linesearch->params = params;

  SLEQP_CALL(sleqp_merit_data_capture(merit_data));

  linesearch->merit_data = merit_data;

  SLEQP_CALL(sleqp_calloc(&linesearch->prod_cache,
                          problem->num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&linesearch->cauchy_point,
                                              problem->num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&linesearch->cauchy_jacobian_prod,
                                              problem->num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&linesearch->newton_jacobian_prod,
                                              problem->num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&linesearch->cauchy_cons_val,
                                              problem->num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&linesearch->combined_cons_val,
                                              problem->num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&linesearch->cauchy_newton_direction,
                                              problem->num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&linesearch->violated_multipliers,
                                              problem->num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&linesearch->test_direction,
                                              problem->num_variables));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_linesearch_set_iterate(SleqpLineSearchData* linesearch,
                                           SleqpIterate* iterate,
                                           double penalty_parameter,
                                           double trust_radius)
{
  linesearch->iterate = iterate;
  linesearch->penalty_parameter = penalty_parameter;
  linesearch->trust_radius = trust_radius;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_linesearch_cauchy_step(SleqpLineSearchData* linesearch,
                                           SleqpSparseVec* direction,
                                           SleqpSparseVec* multipliers,
                                           SleqpSparseVec* hessian_direction,
                                           double* step_length,
                                           double* quadratic_merit_value)
{
  SleqpProblem* problem = linesearch->problem;
  SleqpMeritData* merit_data = linesearch->merit_data;

  const double eps = sleqp_params_get_eps(linesearch->params);
  const double zero_eps = sleqp_params_get_zero_eps(linesearch->params);

  SleqpIterate* iterate = linesearch->iterate;
  const double penalty_parameter = linesearch->penalty_parameter;
  const double trust_radius = linesearch->trust_radius;

  (*quadratic_merit_value) = 0.;

  double exact_merit_value;

  SLEQP_CALL(sleqp_merit_func(merit_data,
                              iterate,
                              penalty_parameter,
                              &exact_merit_value));

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
                                            problem->num_constraints,
                                            zero_eps));
  }

  const double eta = sleqp_params_get_cauchy_eta(linesearch->params);
  const double tau = sleqp_params_get_cauchy_tau(linesearch->params);


  int iteration = 0;

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

      double total_violation;

      SLEQP_CALL(sleqp_violation_one_norm(problem,
                                          linesearch->combined_cons_val,
                                          zero_eps,
                                          &total_violation));

      linear_merit_value += penalty_parameter * total_violation;
    }

    // compute quadratic merit value
    {
      (*quadratic_merit_value) = linear_merit_value + (0.5 * hessian_product);
    }

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

        double func_dual = 1.;

        SLEQP_CALL(sleqp_merit_quadratic(merit_data,
                                         iterate,
                                         &func_dual,
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
    if((exact_merit_value - (*quadratic_merit_value)) >= eta*(exact_merit_value - linear_merit_value))
    {
      break;
    }

    // update products
    SLEQP_CALL(sleqp_sparse_vector_scale(direction, tau));

    hessian_product *= (tau*tau);

    objective_dot *= tau;

    SLEQP_CALL(sleqp_sparse_vector_scale(jacobian_product, tau));

    delta *= tau;
  }

  SLEQP_CALL(sleqp_sparse_vector_scale(hessian_direction, delta));

  assert(iteration != LINESEARCH_MAX_IT);

  sleqp_log_debug("Cauchy line search terminated after %d iterations (step length: %f, quadratic merit: %f)",
                  iteration,
                  delta,
                  (*quadratic_merit_value));

  if(step_length)
  {
    *step_length = delta;
  }

  sleqp_assert_is_leq(sleqp_sparse_vector_norm(direction), trust_radius, eps);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_linesearch_trial_step(SleqpLineSearchData* linesearch,
                                          SleqpSparseVec* cauchy_step,
                                          SleqpSparseVec* cauchy_hessian_step,
                                          const double cauchy_quadratic_merit_value,
                                          SleqpSparseVec* newton_step,
                                          SleqpSparseVec* newton_hessian_step,
                                          SleqpSparseVec* multipliers,
                                          SleqpSparseVec* trial_step,
                                          double* step_length,
                                          double* trial_quadratic_merit_value)
{
  SleqpProblem* problem = linesearch->problem;
  SleqpMeritData* merit_data = linesearch->merit_data;
  SleqpIterate* iterate = linesearch->iterate;
  const double penalty_parameter = linesearch->penalty_parameter;

  const double eps = sleqp_params_get_eps(linesearch->params);
  const double zero_eps = sleqp_params_get_zero_eps(linesearch->params);

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

    SLEQP_CALL(sleqp_func_hess_prod(problem->func,
                                    &func_dual,
                                    cauchy_step,
                                    multipliers,
                                    linesearch->test_direction));

    sleqp_num_assert(sleqp_sparse_vector_eq(cauchy_hessian_step,
                                            linesearch->test_direction,
                                            eps));

    SLEQP_CALL(sleqp_func_hess_prod(problem->func,
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
                                            problem->num_constraints,
                                            zero_eps));

    SLEQP_CALL(sleqp_sparse_matrix_vector_product(sleqp_iterate_get_cons_jac(iterate),
                                                  newton_step,
                                                  linesearch->prod_cache));

    SLEQP_CALL(sleqp_sparse_vector_from_raw(linesearch->newton_jacobian_prod,
                                            linesearch->prod_cache,
                                            problem->num_constraints,
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

  const double eta = sleqp_params_get_linesearch_eta(linesearch->params);
  const double tau = sleqp_params_get_linesearch_tau(linesearch->params);
  const double cutoff_threshold = sleqp_params_get_linesearch_cutoff(linesearch->params);

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
                                     problem->var_lb,
                                     problem->var_ub,
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
                                                     NULL,
                                                     eps));

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

    return SLEQP_OKAY;
  }

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
                                          zero_eps,
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
      double func_dual = 1.;

      double actual_quadratic_merit;

      SLEQP_CALL(sleqp_merit_quadratic(merit_data,
                                       iterate,
                                       &func_dual,
                                       linesearch->test_direction,
                                       multipliers,
                                       penalty_parameter,
                                       &actual_quadratic_merit));

      sleqp_assert_is_eq(quadratic_merit_value,
                         actual_quadratic_merit,
                         eps);
  }

#endif

    // check convergence or abort if the stepsize becomes too small

    if(quadratic_merit_value <= cauchy_quadratic_merit_value + eta*alpha*quadratic_merit_gradient_product)
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
    else if(alpha <= cutoff_threshold)
    {
      (*step_length) = 0.;
      (*trial_quadratic_merit_value) = cauchy_quadratic_merit_value;

      SLEQP_CALL(sleqp_sparse_vector_copy(cauchy_step, trial_step));

      break;
    }

    // Update
    {
      alpha *= tau;
    }

  }

  assert(iteration != LINESEARCH_MAX_IT);

  sleqp_log_debug("Cauchy-Newton line search terminated after %d iterations (step length: %f, quadratic merit: %f)",
                  iteration,
                  alpha,
                  (*trial_quadratic_merit_value));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE linesearch_free(SleqpLineSearchData** star)
{
  SleqpLineSearchData* linesearch = *star;

  SLEQP_CALL(sleqp_sparse_vector_free(&linesearch->test_direction));

  SLEQP_CALL(sleqp_sparse_vector_free(&linesearch->violated_multipliers));
  SLEQP_CALL(sleqp_sparse_vector_free(&linesearch->cauchy_newton_direction));

  SLEQP_CALL(sleqp_sparse_vector_free(&linesearch->combined_cons_val));
  SLEQP_CALL(sleqp_sparse_vector_free(&linesearch->cauchy_cons_val));
  SLEQP_CALL(sleqp_sparse_vector_free(&linesearch->newton_jacobian_prod));
  SLEQP_CALL(sleqp_sparse_vector_free(&linesearch->cauchy_jacobian_prod));
  SLEQP_CALL(sleqp_sparse_vector_free(&linesearch->cauchy_point));
  sleqp_free(&linesearch->prod_cache);

  SLEQP_CALL(sleqp_merit_data_release(&linesearch->merit_data));

  SLEQP_CALL(sleqp_params_release(&linesearch->params));

  sleqp_free(star);

  return SLEQP_OKAY;
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
