#include "newton.h"

#include <math.h>
#include <trlib.h>

#include "fail.h"
#include "cmp.h"
#include "feas.h"
#include "iterate.h"
#include "log.h"
#include "mem.h"
#include "timer.h"
#include "util.h"
#include "working_set.h"

#include "sparse/sparse_matrix.h"

#include "tr/tr_solver.h"
#include "tr/trlib_solver.h"
#include "tr/steihaug_solver.h"

struct SleqpNewtonData
{
  int refcount;
  SleqpProblem* problem;
  SleqpWorkingStep* working_step;

  SleqpParams* params;
  SleqpOptions* options;

  SleqpIterate* iterate;
  SleqpAugJacobian* jacobian;
  double penalty_parameter;

  bool newton_step_in_working_set;

  SleqpSparseVec* gradient;

  SleqpSparseVec* initial_hessian_product;
  SleqpSparseVec* jacobian_product;

  SleqpSparseVec* sparse_cache;
  SleqpSparseVec* tr_step;
  SleqpSparseVec* tr_hessian_product;

  double* dense_cache;

  SleqpTRSolver* tr_solver;

  SleqpTimer* timer;
};

SLEQP_RETCODE sleqp_newton_data_create(SleqpNewtonData** star,
                                       SleqpProblem* problem,
                                       SleqpWorkingStep* step,
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

  data->working_step = step;
  SLEQP_CALL(sleqp_working_step_capture(data->working_step));

  SLEQP_CALL(sleqp_params_capture(params));
  data->params = params;

  SLEQP_CALL(sleqp_options_capture(options));
  data->options = options;

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->gradient,
                                              num_variables));

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
    SLEQP_CALL(sleqp_steihaug_solver_create(&data->tr_solver,
                                            problem,
                                            params,
                                            options));
  }
  else
  {
    assert(tr_solver == SLEQP_TR_SOLVER_TRLIB);

    SLEQP_CALL(sleqp_trlib_solver_create(&data->tr_solver,
                                         problem,
                                         params,
                                         options));
  }

  SLEQP_CALL(sleqp_timer_create(&(data->timer)));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_newton_set_time_limit(SleqpNewtonData* data,
                                          double time_limit)
{
  return sleqp_tr_solver_set_time_limit(data->tr_solver,
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

  SLEQP_CALL(sleqp_working_step_set_iterate(data->working_step,
                                            iterate,
                                            jacobian,
                                            trust_radius));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_newton_add_violated_multipliers(SleqpNewtonData* data,
                                                    SleqpSparseVec* multipliers)
{
  assert(data->iterate);

  SleqpSparseVec* cons_dual = sleqp_iterate_get_cons_dual(data->iterate);

  SleqpSparseVec* violated_cons_mult = sleqp_working_step_get_violated_cons_multipliers(data->working_step);

  const double zero_eps = sleqp_params_get(data->params,
                                           SLEQP_PARAM_ZERO_EPS);

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(cons_dual,
                                            violated_cons_mult,
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

  sleqp_log_debug("Trust region feasibility residuum: %.14e, stationarity residuum: %.14e",
                  SLEQP_MAX(radius_res, projection_res),
                  stationarity_res);

  return SLEQP_OKAY;
}

// compute the EQP gradient. Given as the sum of the
// EQP Hessian with the initial solution, the objective
// function gradient and the violated multipliers
static
SLEQP_RETCODE compute_gradient(SleqpNewtonData* data,
                               SleqpSparseVec* multipliers)
{
  assert(data->iterate);

  SleqpProblem* problem = data->problem;
  SleqpIterate* iterate = data->iterate;

  const double zero_eps = sleqp_params_get(data->params,
                                           SLEQP_PARAM_ZERO_EPS);

  const double penalty_parameter = data->penalty_parameter;

  SleqpSparseMatrix* cons_jac = sleqp_iterate_get_cons_jac(iterate);
  SleqpSparseVec* func_grad = sleqp_iterate_get_func_grad(iterate);

  double one = 1.;

  SleqpSparseVec* initial_step = sleqp_working_step_get_step(data->working_step);

  SleqpSparseVec* violated_cons_mult = sleqp_working_step_get_violated_cons_multipliers(data->working_step);

  SleqpSparseVec* violated_vars_mult = sleqp_working_step_get_violated_vars_multipliers(data->working_step);

  SLEQP_CALL(sleqp_problem_hess_prod(problem,
                                     &one,
                                     initial_step,
                                     multipliers,
                                     data->initial_hessian_product));

  SLEQP_CALL(sleqp_sparse_vector_add(data->initial_hessian_product,
                                     func_grad,
                                     zero_eps,
                                     data->gradient));

  SLEQP_CALL(sleqp_sparse_matrix_trans_vector_product(cons_jac,
                                                      violated_cons_mult,
                                                      zero_eps,
                                                      data->jacobian_product));

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(data->gradient,
                                            data->jacobian_product,
                                            1.,
                                            penalty_parameter,
                                            zero_eps,
                                            data->sparse_cache));

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(data->sparse_cache,
                                            violated_vars_mult,
                                            1.,
                                            penalty_parameter,
                                            zero_eps,
                                            data->gradient));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_newton_compute_step(SleqpNewtonData* data,
                                        SleqpSparseVec* multipliers,
                                        SleqpSparseVec* newton_step)
{
  assert(data->iterate);

  SleqpProblem* problem = data->problem;
  SleqpIterate* iterate = data->iterate;

  SleqpSparseVec* tr_step = data->tr_step;

  SleqpAugJacobian* jacobian = data->jacobian;
  double tr_dual = 0.;

  const double eps = sleqp_params_get(data->params,
                                      SLEQP_PARAM_EPS);

  const double zero_eps = sleqp_params_get(data->params,
                                           SLEQP_PARAM_ZERO_EPS);

  const double reduced_trust_radius = sleqp_working_step_get_reduced_trust_radius(data->working_step);

  SleqpSparseVec* initial_step = sleqp_working_step_get_step(data->working_step);

  // in this case the only feasible solution is the zero vector
  if(sleqp_is_zero(reduced_trust_radius, zero_eps))
  {
    SLEQP_CALL(sleqp_sparse_vector_copy(initial_step, newton_step));

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

  SLEQP_CALL(compute_gradient(data, multipliers));

  SLEQP_CALL(sleqp_tr_solver_solve(data->tr_solver,
                                   jacobian,
                                   multipliers,
                                   data->gradient,
                                   tr_step,
                                   reduced_trust_radius,
                                   &tr_dual));

  SLEQP_CALL(sleqp_sparse_vector_add(tr_step,
                                     initial_step,
                                     zero_eps,
                                     newton_step));

#if !defined(NDEBUG)

  SLEQP_CALL(print_residuals(data,
                             multipliers,
                             data->gradient,
                             tr_step,
                             reduced_trust_radius));

  // Initial direction and trust region direction
  // must be orthogonal
  {
    double direction_dot;

    SLEQP_CALL(sleqp_sparse_vector_dot(tr_step,
                                       initial_step,
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
  SLEQP_CALL(sleqp_tr_solver_current_rayleigh(data->tr_solver,
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

  SLEQP_CALL(sleqp_tr_solver_release(&data->tr_solver));

  sleqp_free(&data->dense_cache);

  SLEQP_CALL(sleqp_sparse_vector_free(&data->tr_hessian_product));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->tr_step));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->sparse_cache));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->jacobian_product));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->initial_hessian_product));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->gradient));

  SLEQP_CALL(sleqp_aug_jacobian_release(&data->jacobian));
  SLEQP_CALL(sleqp_iterate_release(&data->iterate));

  SLEQP_CALL(sleqp_options_release(&data->options));
  SLEQP_CALL(sleqp_params_release(&data->params));

  SLEQP_CALL(sleqp_working_step_release(&data->working_step));

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
