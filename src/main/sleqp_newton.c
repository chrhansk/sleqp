#include "sleqp_newton.h"

#include <math.h>
#include <trlib.h>

#include "sleqp_assert.h"
#include "sleqp_cmp.h"
#include "sleqp_iterate.h"
#include "sleqp_log.h"
#include "sleqp_mem.h"
#include "sleqp_timer.h"
#include "sleqp_util.h"
#include "sleqp_tr_solver.h"

struct SleqpNewtonData
{
  int refcount;
  SleqpProblem* problem;
  SleqpParams* params;
  SleqpOptions* options;

  SleqpSparseVec* initial_solution;
  SleqpSparseVec* initial_cons_val;

  SleqpSparseVec* lower_diff;
  SleqpSparseVec* upper_diff;

  SleqpSparseVec* violated_multipliers;
  SleqpSparseVec* multipliers;

  SleqpSparseVec* gradient;

  SleqpSparseVec* initial_rhs;
  SleqpSparseVec* initial_hessian_product;
  SleqpSparseVec* jacobian_product;

  SleqpSparseVec* sparse_cache;
  double* dense_cache;

  SleqpTRSolver* trust_region_solver;
};

SLEQP_RETCODE sleqp_newton_data_create(SleqpNewtonData** star,
                                       SleqpProblem* problem,
                                       SleqpParams* params,
                                       SleqpOptions* options)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpNewtonData* data = *star;

  *data = (SleqpNewtonData) {0};

  data->refcount = 1;

  data->problem = problem;
  data->params = params;
  data->options = options;

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->initial_solution,
                                              problem->num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->initial_cons_val,
                                              problem->num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->lower_diff,
                                              problem->num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->upper_diff,
                                              problem->num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->violated_multipliers,
                                              problem->num_constraints));
  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->multipliers,
                                              problem->num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->gradient,
                                              problem->num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->initial_rhs,
                                              problem->num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->initial_hessian_product,
                                              problem->num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->jacobian_product,
                                              problem->num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->sparse_cache,
                                              problem->num_variables));

  SLEQP_CALL(sleqp_calloc(&data->dense_cache,
                          SLEQP_MAX(problem->num_variables, problem->num_constraints)));

  SLEQP_CALL(sleqp_tr_solver_create(&data->trust_region_solver,
                                    problem,
                                    params,
                                    options));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE get_initial_rhs(SleqpNewtonData* data,
                                     SleqpIterate* iterate,
                                     SleqpAugJacobian* jacobian)
{
  SleqpProblem* problem = data->problem;

  SleqpSparseVec* initial_rhs = data->initial_rhs;
  SleqpWorkingSet* working_set = sleqp_iterate_get_working_set(iterate);

  const double eps = sleqp_params_get_eps(data->params);

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
    SleqpSparseVec* var_lb = problem->var_lb;
    SleqpSparseVec* var_ub = problem->var_ub;

    SleqpSparseVec* lower_diff = data->lower_diff;
    SleqpSparseVec* upper_diff = data->upper_diff;

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(values, var_ub, -1., 1., eps, upper_diff));

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(values, var_lb, -1., 1., eps, lower_diff));

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
    SleqpSparseVec* cons_lb = problem->cons_lb;
    SleqpSparseVec* cons_ub = problem->cons_ub;

    SleqpSparseVec* lower_diff = data->lower_diff;
    SleqpSparseVec* upper_diff = data->upper_diff;

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(values, cons_ub, -1., 1., eps, upper_diff));

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(values, cons_lb, -1., 1., eps, lower_diff));

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

static SLEQP_RETCODE matrix_pop_column(SleqpSparseMatrix* matrix)
{
  const int num_cols = sleqp_sparse_matrix_get_num_cols(matrix);

  assert(num_cols > 0);

  const int last_col = num_cols - 1;

  SLEQP_CALL(sleqp_sparse_matrix_pop_column(matrix, last_col));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_newton_set_time_limit(SleqpNewtonData* data,
                                          double time_limit)
{
  return sleqp_tr_solver_set_time_limit(data->trust_region_solver,
                                        time_limit);
}

SLEQP_RETCODE sleqp_newton_compute_step(SleqpNewtonData* data,
                                        SleqpIterate* iterate,
                                        SleqpAugJacobian* jacobian,
                                        double trust_radius,
                                        double penalty_parameter,
                                        SleqpSparseVec* newton_step)
{
  SleqpProblem* problem = data->problem;

  const double eps = sleqp_params_get_eps(data->params);
  const double zero_eps = sleqp_params_get_zero_eps(data->params);

  SLEQP_CALL(get_initial_rhs(data, iterate, jacobian));

  SLEQP_CALL(sleqp_aug_jacobian_min_norm_solution(jacobian,
                                                  data->initial_rhs,
                                                  data->initial_solution));

#if !defined(NDEBUG)

  // Initial direction must be in working set
  {
    bool in_working_set = false;

    SLEQP_CALL(sleqp_direction_in_working_set(problem,
                                              iterate,
                                              data->initial_solution,
                                              data->dense_cache,
                                              eps,
                                              &in_working_set));

    sleqp_num_assert_msg(in_working_set,
                         "Initial Newton direction must respect working set linearization");
  }

#endif

  const double norm_ratio = .8;
  bool newton_step_in_working_set = true;

  // rescale min norm solution if required.
  {
    double unscaled_norm = sleqp_sparse_vector_norm(data->initial_solution);

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
        newton_step_in_working_set = false;

        SLEQP_CALL(sleqp_sparse_vector_scale(data->initial_solution, alpha));

        // we know that the scaled initial solution
        // has norm equal to norm_ratio * trust_radius
        trust_radius *= sqrt(1. - norm_ratio * norm_ratio);
      }
    }
  }

  SleqpSparseMatrix* cons_jac = sleqp_iterate_get_cons_jac(iterate);
  SleqpSparseVec* func_grad = sleqp_iterate_get_func_grad(iterate);

  SleqpWorkingSet* working_set = sleqp_iterate_get_working_set(iterate);
  int working_set_size = sleqp_working_set_size(working_set);

  // in either case the only feasible solution is the zero vector
  if(sleqp_is_zero(trust_radius, zero_eps) || problem->num_variables == working_set_size)
  {
    SLEQP_CALL(sleqp_sparse_vector_copy(data->initial_solution, newton_step));

    return SLEQP_OKAY;
  }

  // compute the EQP multipliers. Given as the sum of the
  // estimated multipliers and penalties on violated
  // constraints
  {
    // First: Compute linearized constraint values according to initial solution

    SLEQP_CALL(sleqp_sparse_matrix_vector_product(cons_jac,
                                                  data->initial_solution,
                                                  data->dense_cache));

    SLEQP_CALL(sleqp_sparse_vector_from_raw(data->sparse_cache,
                                            data->dense_cache,
                                            problem->num_constraints,
                                            zero_eps));

    SLEQP_CALL(sleqp_sparse_vector_add(sleqp_iterate_get_cons_val(iterate),
                                       data->sparse_cache,
                                       eps,
                                       data->initial_cons_val));

    // Second: Compute multipliers for violated constraints
    SLEQP_CALL(sleqp_get_violated_multipliers(problem,
                                              data->initial_cons_val,
                                              data->violated_multipliers,
                                              working_set,
                                              zero_eps));

    // Finally: Add the iterate multipliers to the given ones
    SLEQP_CALL(sleqp_sparse_vector_add_scaled(sleqp_iterate_get_cons_dual(iterate),
                                              data->violated_multipliers,
                                              1.,
                                              penalty_parameter,
                                              zero_eps,
                                              data->multipliers));
  }

  // compute the EQP gradient. Given as the sum of the
  // EQP Hessian with the initial solution, the objective
  // function gradient and the violated multipliers
  {
    double one = 1.;

    SLEQP_CALL(sleqp_func_hess_prod(problem->func,
                                    &one,
                                    data->initial_solution,
                                    data->multipliers,
                                    data->initial_hessian_product));

    SLEQP_CALL(sleqp_sparse_vector_add(data->initial_hessian_product,
                                       func_grad,
                                       eps,
                                       data->sparse_cache));

    SLEQP_CALL(sleqp_sparse_matrix_trans_vector_product(cons_jac,
                                                        data->multipliers,
                                                        eps,
                                                        data->jacobian_product));

    SLEQP_CALL(sleqp_sparse_vector_add(data->sparse_cache,
                                       data->jacobian_product,
                                       eps,
                                       data->gradient));
  }

  SLEQP_CALL(sleqp_tr_solver_solve(data->trust_region_solver,
                                   jacobian,
                                   data->multipliers,
                                   data->gradient,
                                   data->sparse_cache,
                                   trust_radius));

  SLEQP_CALL(sleqp_sparse_vector_add(data->sparse_cache,
                                     data->initial_solution,
                                     zero_eps,
                                     newton_step));

#if !defined(NDEBUG)

  // Initial direction and trust region direction
  // must be orthogonal
  {
    double direction_dot;
    sleqp_sparse_vector_dot(data->sparse_cache,
                            data->initial_solution,
                            &direction_dot);

    sleqp_assert_is_zero(direction_dot, eps);
  }

  if(newton_step_in_working_set)
  {
    // Direction must be in working set
    bool in_working_set = false;

    SLEQP_CALL(sleqp_direction_in_working_set(problem,
                                              iterate,
                                              newton_step,
                                              data->dense_cache,
                                              eps,
                                              &in_working_set));

    sleqp_num_assert_msg(in_working_set,
                         "Newton step must respect working set linearization");
  }

#endif

  return SLEQP_OKAY;
}

static SLEQP_RETCODE newton_data_free(SleqpNewtonData** star)
{
  SleqpNewtonData* data = *star;

  if(!data)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_tr_solver_release(&data->trust_region_solver));

  sleqp_free(&data->dense_cache);

  SLEQP_CALL(sleqp_sparse_vector_free(&data->sparse_cache));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->jacobian_product));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->initial_hessian_product));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->initial_rhs));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->gradient));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->multipliers));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->violated_multipliers));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->upper_diff));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->lower_diff));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->initial_cons_val));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->initial_solution));

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
