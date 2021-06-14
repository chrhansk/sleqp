#include "sleqp_newton.h"

#include <math.h>
#include <trlib.h>

#include "sleqp_assert.h"
#include "sleqp_cmp.h"
#include "sleqp_feas.h"
#include "sleqp_iterate.h"
#include "sleqp_log.h"
#include "sleqp_mem.h"
#include "sleqp_timer.h"
#include "sleqp_util.h"

#include "tr/sleqp_tr_solver.h"
#include "tr/sleqp_trlib_solver.h"
#include "tr/sleqp_steihaug_solver.h"

struct SleqpNewtonData
{
  int refcount;
  SleqpProblem* problem;
  SleqpParams* params;
  SleqpOptions* options;

  SleqpIterate* iterate;
  SleqpAugJacobian* jacobian;
  double trust_radius;
  double penalty_parameter;

  SleqpSparseVec* initial_direction;
  SleqpSparseVec* initial_point;

  SleqpSparseVec* initial_cons_val;

  SleqpSparseVec* lower_diff;
  SleqpSparseVec* upper_diff;

  SleqpSparseVec* violated_variable_multipliers;
  SleqpSparseVec* violated_constraint_multipliers;

  bool newton_step_in_working_set;

  SleqpSparseVec* gradient;

  SleqpSparseVec* initial_rhs;
  SleqpSparseVec* initial_hessian_product;
  SleqpSparseVec* jacobian_product;

  SleqpSparseVec* sparse_cache;
  SleqpSparseVec* tr_step;
  SleqpSparseVec* tr_hessian_product;

  double* dense_cache;

  SleqpTRSolver* trust_region_solver;

  SleqpTimer* timer;
};

SLEQP_RETCODE sleqp_newton_data_create(SleqpNewtonData** star,
                                       SleqpProblem* problem,
                                       SleqpParams* params,
                                       SleqpOptions* options)
{
  SLEQP_CALL(sleqp_malloc(star));

  const int num_variables = sleqp_problem_num_variables(problem);
  const int num_constraints = sleqp_problem_num_constraints(problem);

  SleqpNewtonData* data = *star;

  *data = (SleqpNewtonData) {0};

  data->refcount = 1;

  data->problem = problem;
  SLEQP_CALL(sleqp_problem_capture(data->problem));

  SLEQP_CALL(sleqp_params_capture(params));
  data->params = params;

  SLEQP_CALL(sleqp_options_capture(options));
  data->options = options;

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->initial_direction,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->initial_point,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->initial_cons_val,
                                              num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->lower_diff,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->upper_diff,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->violated_variable_multipliers,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->violated_constraint_multipliers,
                                              num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->gradient,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->initial_rhs,
                                              num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->initial_hessian_product,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->jacobian_product,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->sparse_cache,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->tr_step,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->tr_hessian_product,
                                              num_variables));

  SLEQP_CALL(sleqp_alloc_array(&data->dense_cache,
                               SLEQP_MAX(num_variables, num_constraints)));

  SLEQP_TR_SOLVER tr_solver = sleqp_options_get_int(options, SLEQP_OPTION_INT_TR_SOLVER);

  if(tr_solver == SLEQP_TR_SOLVER_AUTO)
  {
    SleqpFunc* func = sleqp_problem_func(problem);

    if(sleqp_func_has_psd_hessian(func))
    {
      tr_solver = SLEQP_TR_SOLVER_CG;
    }
    else
    {
      tr_solver = SLEQP_TR_SOLVER_TRLIB;
    }
  }

  if(tr_solver == SLEQP_TR_SOLVER_CG)
  {
    SLEQP_CALL(sleqp_steihaug_solver_create(&data->trust_region_solver,
                                            problem,
                                            params,
                                            options));
  }
  else
  {
    assert(tr_solver == SLEQP_TR_SOLVER_TRLIB);

    SLEQP_CALL(sleqp_trlib_solver_create(&data->trust_region_solver,
                                         problem,
                                         params,
                                         options));
  }

  SLEQP_CALL(sleqp_timer_create(&(data->timer)));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE get_initial_rhs(SleqpNewtonData* data,
                                     SleqpIterate* iterate,
                                     SleqpAugJacobian* jacobian)
{
  SleqpProblem* problem = data->problem;

  SleqpSparseVec* initial_rhs = data->initial_rhs;
  SleqpWorkingSet* working_set = sleqp_iterate_get_working_set(iterate);

  const double eps = sleqp_params_get(data->params, SLEQP_PARAM_EPS);

  const double zero_eps = sleqp_params_get(data->params, SLEQP_PARAM_ZERO_EPS);

  const int working_set_size = sleqp_working_set_size(working_set);

  SLEQP_CALL(sleqp_sparse_vector_resize(initial_rhs, working_set_size));

  {
    SLEQP_CALL(sleqp_sparse_vector_clear(initial_rhs));

    SLEQP_CALL(sleqp_sparse_vector_reserve(initial_rhs, working_set_size));

    SLEQP_CALL(sleqp_sparse_vector_resize(initial_rhs, working_set_size));
  }

  // variables
  {
    SleqpSparseVec* values = sleqp_iterate_get_primal(iterate);
    SleqpSparseVec* var_lb = sleqp_problem_var_lb(problem);
    SleqpSparseVec* var_ub = sleqp_problem_var_ub(problem);

    SleqpSparseVec* lower_diff = data->lower_diff;
    SleqpSparseVec* upper_diff = data->upper_diff;

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(values, var_ub, -1., 1., zero_eps, upper_diff));

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(values, var_lb, -1., 1., zero_eps, lower_diff));

    int k_lower = 0, k_upper = 0;

    while(k_lower < lower_diff->nnz || k_upper < upper_diff->nnz)
    {
      bool valid_lower = (k_lower < lower_diff->nnz);
      bool valid_upper = (k_upper < upper_diff->nnz);

      const int i_lower = valid_lower ? lower_diff->indices[k_lower] : lower_diff->dim + 1;
      const int i_upper = valid_upper ? upper_diff->indices[k_upper] : upper_diff->dim + 1;

      const int i_combined = SLEQP_MIN(i_lower, i_upper);

      valid_lower = valid_lower && (i_lower == i_combined);
      valid_upper = valid_upper && (i_upper == i_combined);

      const double lower_value = valid_lower ? lower_diff->data[k_lower] : 0.;
      const double upper_value = valid_upper ? upper_diff->data[k_upper] : 0.;

      const int i_set = sleqp_working_set_get_variable_index(working_set, i_combined);

      const SLEQP_ACTIVE_STATE var_state = sleqp_working_set_get_variable_state(working_set,
                                                                                i_combined);

      assert(var_state == SLEQP_INACTIVE || i_set != SLEQP_NONE);

      if(var_state == SLEQP_ACTIVE_UPPER)
      {
        SLEQP_CALL(sleqp_sparse_vector_push(initial_rhs,
                                            i_set,
                                            upper_value));
      }
      else if(var_state == SLEQP_ACTIVE_LOWER)
      {
        SLEQP_CALL(sleqp_sparse_vector_push(initial_rhs,
                                            i_set,
                                            lower_value));
      }
      else if(var_state == SLEQP_ACTIVE_BOTH)
      {
        sleqp_assert_is_eq(lower_value, upper_value, eps);

        SLEQP_CALL(sleqp_sparse_vector_push(initial_rhs,
                                            i_set,
                                            lower_value));
      }

      if(i_lower == i_combined)
      {
        ++k_lower;
      }

      if(i_upper == i_combined)
      {
        ++k_upper;
      }

    }
  }

  // constraints
  {
    SleqpSparseVec* values = sleqp_iterate_get_cons_val(iterate);
    SleqpSparseVec* cons_lb = sleqp_problem_cons_lb(problem);
    SleqpSparseVec* cons_ub = sleqp_problem_cons_ub(problem);

    SleqpSparseVec* lower_diff = data->lower_diff;
    SleqpSparseVec* upper_diff = data->upper_diff;

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(values, cons_ub, -1., 1., zero_eps, upper_diff));

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(values, cons_lb, -1., 1., zero_eps, lower_diff));

    int k_lower = 0, k_upper = 0;

    while(k_lower < lower_diff->nnz || k_upper < upper_diff->nnz)
    {
      bool valid_lower = (k_lower < lower_diff->nnz);
      bool valid_upper = (k_upper < upper_diff->nnz);

      const int i_lower = valid_lower ? lower_diff->indices[k_lower] : lower_diff->dim + 1;
      const int i_upper = valid_upper ? upper_diff->indices[k_upper] : upper_diff->dim + 1;

      const int i_combined = SLEQP_MIN(i_lower, i_upper);

      valid_lower = valid_lower && (i_lower == i_combined);
      valid_upper = valid_upper && (i_upper == i_combined);

      const double lower_value = valid_lower ? lower_diff->data[k_lower] : 0.;
      const double upper_value = valid_upper ? upper_diff->data[k_upper] : 0.;

      const int i_set = sleqp_working_set_get_constraint_index(working_set, i_combined);

      const SLEQP_ACTIVE_STATE cons_state = sleqp_working_set_get_constraint_state(working_set,
                                                                                   i_combined);

      assert(cons_state == SLEQP_INACTIVE || i_set != SLEQP_NONE);

      if(cons_state == SLEQP_ACTIVE_UPPER)
      {
        SLEQP_CALL(sleqp_sparse_vector_push(initial_rhs,
                                            i_set,
                                            upper_value));
      }
      else if(cons_state == SLEQP_ACTIVE_LOWER)
      {
        SLEQP_CALL(sleqp_sparse_vector_push(initial_rhs,
                                            i_set,
                                            lower_value));
      }
      else if(cons_state == SLEQP_ACTIVE_BOTH)
      {
        sleqp_assert_is_eq(lower_value, upper_value, eps);

        SLEQP_CALL(sleqp_sparse_vector_push(initial_rhs,
                                            i_set,
                                            lower_value));
      }

      if(valid_lower)
      {
        ++k_lower;
      }

      if(valid_upper)
      {
        ++k_upper;
      }
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_newton_set_time_limit(SleqpNewtonData* data,
                                          double time_limit)
{
  return sleqp_tr_solver_set_time_limit(data->trust_region_solver,
                                        time_limit);
}

SleqpTimer* sleqp_newton_get_timer(SleqpNewtonData* data)
{
  return data->timer;
}

SLEQP_RETCODE sleqp_newton_set_iterate(SleqpNewtonData* data,
                                       SleqpIterate* iterate,
                                       SleqpAugJacobian* jacobian,
                                       double trust_radius,
                                       double penalty_parameter)
{
  data->trust_radius = trust_radius;
  data->penalty_parameter = penalty_parameter;

  {
    SLEQP_CALL(sleqp_iterate_release(&data->iterate));

    SLEQP_CALL(sleqp_iterate_capture(iterate));

    data->iterate = iterate;
  }

  {
    SLEQP_CALL(sleqp_aug_jacobian_release(&data->jacobian));

    SLEQP_CALL(sleqp_aug_jacobian_capture(jacobian));

    data->jacobian = jacobian;
  }

  SleqpProblem* problem = data->problem;

  const int num_constraints = sleqp_problem_num_constraints(problem);

  const double eps = sleqp_params_get(data->params,
                                      SLEQP_PARAM_EPS);

  const double zero_eps = sleqp_params_get(data->params,
                                           SLEQP_PARAM_ZERO_EPS);

  SLEQP_CALL(get_initial_rhs(data, iterate, jacobian));

  SLEQP_CALL(sleqp_aug_jacobian_min_norm_solution(jacobian,
                                                  data->initial_rhs,
                                                  data->initial_direction));

#if !defined(NDEBUG)

  // Initial direction must be in working set
  {
    bool in_working_set = false;

    SLEQP_CALL(sleqp_direction_in_working_set(problem,
                                              iterate,
                                              data->initial_direction,
                                              data->dense_cache,
                                              eps,
                                              &in_working_set));

    sleqp_num_assert(in_working_set);
  }

#endif

  const double norm_ratio = .8;

  data->newton_step_in_working_set = true;

  // rescale min norm solution if required.
  {
    double unscaled_norm = sleqp_sparse_vector_norm(data->initial_direction);

    double initial_norm = unscaled_norm;

    double alpha = 1.;

    if(initial_norm != 0.)
    {
      assert(initial_norm > 0.);

      alpha = (norm_ratio * trust_radius) / (initial_norm);

      alpha = SLEQP_MIN(alpha, 1.);

      if(sleqp_is_eq(alpha, 1., eps))
      {
        // no scaling required...

        const double initial_norm_sq = initial_norm * initial_norm;

        const double trust_radius_sq = trust_radius * trust_radius;

        //sleqp_assert_is_lt(initial_norm_sq, trust_radius_sq, eps);

        trust_radius = sqrt(trust_radius_sq - initial_norm_sq);
      }
      else
      {
        data->newton_step_in_working_set = false;

        SLEQP_CALL(sleqp_sparse_vector_scale(data->initial_direction, alpha));

        // we know that the scaled initial solution
        // has norm equal to norm_ratio * trust_radius
        data->trust_radius *= sqrt(1. - norm_ratio * norm_ratio);
      }
    }
  }

  SLEQP_CALL(sleqp_sparse_vector_add(sleqp_iterate_get_primal(iterate),
                                     data->initial_direction,
                                     zero_eps,
                                     data->initial_point));

  SleqpSparseMatrix* cons_jac = sleqp_iterate_get_cons_jac(iterate);
  SleqpWorkingSet* working_set = sleqp_iterate_get_working_set(iterate);

  // Compute linearized constraint values at initial direction
  {
    SLEQP_CALL(sleqp_sparse_matrix_vector_product(cons_jac,
                                                  data->initial_direction,
                                                  data->dense_cache));

    SLEQP_CALL(sleqp_sparse_vector_from_raw(data->sparse_cache,
                                            data->dense_cache,
                                            num_constraints,
                                            zero_eps));

    SLEQP_CALL(sleqp_sparse_vector_add(sleqp_iterate_get_cons_val(iterate),
                                       data->sparse_cache,
                                       zero_eps,
                                       data->initial_cons_val));
  }

  // Compute violated multipliers
  {
    SLEQP_CALL(sleqp_violated_variable_multipliers(problem,
                                                   data->initial_point,
                                                   data->violated_variable_multipliers,
                                                   working_set));

    SLEQP_CALL(sleqp_violated_constraint_multipliers(problem,
                                                     data->initial_cons_val,
                                                     data->violated_constraint_multipliers,
                                                     working_set));
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_newton_compute_multipliers(SleqpNewtonData* data,
                                               SleqpSparseVec* multipliers)
{
  assert(data->iterate);

  SleqpSparseVec* cons_dual = sleqp_iterate_get_cons_dual(data->iterate);

  const double zero_eps = sleqp_params_get(data->params,
                                           SLEQP_PARAM_ZERO_EPS);

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(cons_dual,
                                            data->violated_constraint_multipliers,
                                            1.,
                                            data->penalty_parameter,
                                            zero_eps,
                                            multipliers));

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE print_residuals(SleqpNewtonData* data,
                              const SleqpSparseVec* multipliers,
                              const SleqpSparseVec* gradient,
                              SleqpSparseVec* tr_step,
                              double trust_radius)
{
  SleqpProblem* problem = data->problem;

  SleqpAugJacobian* jacobian = data->jacobian;
  SleqpSparseVec* sparse_cache = data->sparse_cache;
  SleqpSparseVec* tr_prod = data->tr_hessian_product;

  const double zero_eps = sleqp_params_get(data->params,
                                           SLEQP_PARAM_ZERO_EPS);

  const double step_norm = sleqp_sparse_vector_norm(tr_step);
  const double radius_res = step_norm - trust_radius;

  SLEQP_CALL(sleqp_aug_jacobian_projection(jacobian,
                                           tr_step,
                                           sparse_cache,
                                           NULL));

  const double projection_res = sleqp_sparse_vector_inf_norm(sparse_cache);

  double one = 1.;

  SLEQP_CALL(sleqp_problem_hess_prod(problem,
                                     &one,
                                     tr_step,
                                     multipliers,
                                     tr_prod));

  SLEQP_CALL(sleqp_sparse_vector_add(tr_prod,
                                     gradient,
                                     zero_eps,
                                     sparse_cache));

  SLEQP_CALL(sleqp_aug_jacobian_projection(jacobian,
                                           sparse_cache,
                                           tr_prod,
                                           NULL));

  const double stationarity_res = sleqp_sparse_vector_inf_norm(tr_prod);

  sleqp_log_debug("Trust region feasibility residuum: %f, stationarity residuum: %f",
                  SLEQP_MAX(radius_res, projection_res),
                  stationarity_res);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_newton_compute_step(SleqpNewtonData* data,
                                        SleqpSparseVec* multipliers,
                                        SleqpSparseVec* newton_step)
{
  assert(data->iterate);

  SleqpProblem* problem = data->problem;
  SleqpIterate* iterate = data->iterate;

  SleqpSparseMatrix* cons_jac = sleqp_iterate_get_cons_jac(iterate);
  SleqpSparseVec* func_grad = sleqp_iterate_get_func_grad(iterate);
  SleqpSparseVec* tr_step = data->tr_step;

  SleqpAugJacobian* jacobian = data->jacobian;
  double tr_dual = 0.;

  const double eps = sleqp_params_get(data->params,
                                      SLEQP_PARAM_EPS);

  const double zero_eps = sleqp_params_get(data->params,
                                           SLEQP_PARAM_ZERO_EPS);

  const double trust_radius = data->trust_radius;
  const double penalty_parameter = data->penalty_parameter;

  // in this case the only feasible solution is the zero vector
  if(sleqp_is_zero(trust_radius, zero_eps))
  {
    SLEQP_CALL(sleqp_sparse_vector_copy(data->initial_direction, newton_step));

    return SLEQP_OKAY;
  }

  /*
    SleqpWorkingSet* working_set = sleqp_iterate_get_working_set(iterate);
    int working_set_size = sleqp_working_set_size(working_set);
  */

  SleqpFunc* func = sleqp_problem_func(problem);

  SleqpTimer* hess_timer = sleqp_func_get_hess_timer(func);

  SleqpTimer* subst_timer = sleqp_aug_jacobian_get_substitution_timer(data->jacobian);

  const double hess_before = sleqp_timer_elapsed(hess_timer);
  const double subst_before = sleqp_timer_elapsed(subst_timer);

  SLEQP_CALL(sleqp_timer_start(data->timer));

  // compute the EQP gradient. Given as the sum of the
  // EQP Hessian with the initial solution, the objective
  // function gradient and the violated multipliers
  {
    double one = 1.;

    SLEQP_CALL(sleqp_problem_hess_prod(problem,
                                       &one,
                                       data->initial_direction,
                                       multipliers,
                                       data->initial_hessian_product));

    SLEQP_CALL(sleqp_sparse_vector_add(data->initial_hessian_product,
                                       func_grad,
                                       zero_eps,
                                       data->gradient));

    SLEQP_CALL(sleqp_sparse_matrix_trans_vector_product(cons_jac,
                                                        data->violated_constraint_multipliers,
                                                        zero_eps,
                                                        data->jacobian_product));

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(data->gradient,
                                              data->jacobian_product,
                                              1.,
                                              penalty_parameter,
                                              zero_eps,
                                              data->sparse_cache));

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(data->sparse_cache,
                                              data->violated_variable_multipliers,
                                              1.,
                                              penalty_parameter,
                                              zero_eps,
                                              data->gradient));
  }

  SLEQP_CALL(sleqp_tr_solver_solve(data->trust_region_solver,
                                   jacobian,
                                   multipliers,
                                   data->gradient,
                                   tr_step,
                                   trust_radius,
                                   &tr_dual));

  SLEQP_CALL(sleqp_sparse_vector_add(tr_step,
                                     data->initial_direction,
                                     zero_eps,
                                     newton_step));

#if !defined(NDEBUG)

  SLEQP_CALL(print_residuals(data,
                             multipliers,
                             data->gradient,
                             tr_step,
                             trust_radius));

  // Initial direction and trust region direction
  // must be orthogonal
  {
    double direction_dot;

    SLEQP_CALL(sleqp_sparse_vector_dot(tr_step,
                                       data->initial_direction,
                                       &direction_dot));

    sleqp_assert_is_zero(direction_dot, eps);
  }

  if(data->newton_step_in_working_set)
  {
    // Direction must be in working set
    bool in_working_set = false;

    SLEQP_CALL(sleqp_direction_in_working_set(problem,
                                              iterate,
                                              newton_step,
                                              data->dense_cache,
                                              eps,
                                              &in_working_set));

    sleqp_num_assert(in_working_set);
  }

#endif

  SLEQP_CALL(sleqp_timer_stop(data->timer));

  const double hess_elapsed = sleqp_timer_elapsed(hess_timer) - hess_before;
  const double subst_elapsed = sleqp_timer_elapsed(subst_timer) - subst_before;

  SLEQP_CALL(sleqp_timer_add(data->timer, -(hess_elapsed + subst_elapsed)));

  return SLEQP_OKAY;
}


SLEQP_RETCODE sleqp_newton_current_rayleigh(SleqpNewtonData* data,
                                            double* min_rayleigh,
                                            double* max_rayleigh)
{
  SLEQP_CALL(sleqp_tr_solver_current_rayleigh(data->trust_region_solver,
                                              min_rayleigh,
                                              max_rayleigh));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE newton_data_free(SleqpNewtonData** star)
{
  SleqpNewtonData* data = *star;

  if(!data)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_timer_free(&(data->timer)));

  SLEQP_CALL(sleqp_tr_solver_release(&data->trust_region_solver));

  sleqp_free(&data->dense_cache);

  SLEQP_CALL(sleqp_sparse_vector_free(&data->tr_hessian_product));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->tr_step));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->sparse_cache));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->jacobian_product));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->initial_hessian_product));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->initial_rhs));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->gradient));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->violated_constraint_multipliers));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->violated_variable_multipliers));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->upper_diff));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->lower_diff));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->initial_cons_val));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->initial_point));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->initial_direction));

  SLEQP_CALL(sleqp_aug_jacobian_release(&data->jacobian));
  SLEQP_CALL(sleqp_iterate_release(&data->iterate));

  SLEQP_CALL(sleqp_options_release(&data->options));
  SLEQP_CALL(sleqp_params_release(&data->params));

  SLEQP_CALL(sleqp_problem_release(&data->problem));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_newton_data_capture(SleqpNewtonData* data)
{
  ++data->refcount;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_newton_data_release(SleqpNewtonData** star)
{
  SleqpNewtonData* data = *star;

  if(!data)
  {
    return SLEQP_OKAY;
  }

  if(--data->refcount == 0)
  {
    SLEQP_CALL(newton_data_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
