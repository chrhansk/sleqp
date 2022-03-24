#include "merit.h"

#include "cmp.h"
#include "feas.h"
#include "log.h"
#include "mem.h"
#include "util.h"

struct SleqpMeritData
{
  int refcount;

  SleqpProblem* problem;
  SleqpParams* params;

  double* dense_cache;
  int cache_size;

  SleqpVec* jac_dot_sparse;
  SleqpVec* combined_cons_val;

  SleqpVec* multipliers;
  SleqpVec* sparse_cache;
};

SLEQP_RETCODE
sleqp_merit_create(SleqpMerit** star,
                   SleqpProblem* problem,
                   SleqpParams* params)
{
  SLEQP_CALL(sleqp_malloc(star));

  const int num_variables   = sleqp_problem_num_vars(problem);
  const int num_constraints = sleqp_problem_num_cons(problem);

  SleqpMerit* merit = *star;

  merit->refcount = 1;

  merit->problem = problem;
  SLEQP_CALL(sleqp_problem_capture(merit->problem));

  SLEQP_CALL(sleqp_params_capture(params));
  merit->params = params;

  merit->cache_size = SLEQP_MAX(num_constraints, num_variables);

  SLEQP_CALL(sleqp_alloc_array(&merit->dense_cache, merit->cache_size));

  SLEQP_CALL(sleqp_vec_create_empty(&merit->jac_dot_sparse, num_constraints));

  SLEQP_CALL(
    sleqp_vec_create_empty(&merit->combined_cons_val, num_constraints));

  SLEQP_CALL(sleqp_vec_create_empty(&merit->multipliers, num_constraints));

  SLEQP_CALL(sleqp_vec_create_empty(&merit->sparse_cache, num_constraints));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_merit_func(SleqpMerit* merit,
                 SleqpIterate* iterate,
                 double penalty_parameter,
                 double* merit_value)
{
  SleqpProblem* problem = merit->problem;

  *merit_value = sleqp_iterate_obj_val(iterate);

  double total_violation;

  SLEQP_CALL(sleqp_violation_one_norm(problem,
                                      sleqp_iterate_cons_val(iterate),
                                      &total_violation));

  (*merit_value) += penalty_parameter * total_violation;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_merit_linear(SleqpMerit* merit,
                   SleqpIterate* iterate,
                   const SleqpVec* direction,
                   double penalty_parameter,
                   double* merit_value)
{
  SleqpProblem* problem = merit->problem;

  const int num_constraints = sleqp_problem_num_cons(problem);

  const double zero_eps
    = sleqp_params_value(merit->params, SLEQP_PARAM_ZERO_EPS);

  double objective_dot;

  SLEQP_CALL(
    sleqp_vec_dot(sleqp_iterate_obj_grad(iterate), direction, &objective_dot));

  (*merit_value) = sleqp_iterate_obj_val(iterate) + objective_dot;

  SLEQP_CALL(sleqp_sparse_matrix_vector_product(sleqp_iterate_cons_jac(iterate),
                                                direction,
                                                merit->dense_cache));

  SLEQP_CALL(sleqp_vec_set_from_raw(merit->jac_dot_sparse,
                                    merit->dense_cache,
                                    num_constraints,
                                    zero_eps));

  SLEQP_CALL(sleqp_vec_add(merit->jac_dot_sparse,
                           sleqp_iterate_cons_val(iterate),
                           zero_eps,
                           merit->combined_cons_val));

  double total_violation;

  SLEQP_CALL(sleqp_violation_one_norm(problem,
                                      merit->combined_cons_val,
                                      &total_violation));

  (*merit_value) += penalty_parameter * total_violation;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_merit_quadratic(SleqpMerit* merit,
                      SleqpIterate* iterate,
                      const double* obj_dual,
                      const SleqpVec* direction,
                      const SleqpVec* cons_duals,
                      double penalty_parameter,
                      double* merit_value)
{
  SleqpProblem* problem = merit->problem;

  double linear_merit_value;

  SLEQP_CALL(sleqp_merit_linear(merit,
                                iterate,
                                direction,
                                penalty_parameter,
                                &linear_merit_value));

  double bilinear_product;

  SLEQP_CALL(sleqp_problem_hess_bilinear(problem,
                                         obj_dual,
                                         direction,
                                         cons_duals,
                                         &bilinear_product));

  (*merit_value) = linear_merit_value + (0.5 * bilinear_product);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
merit_free(SleqpMerit** star)
{
  SleqpMerit* merit = *star;

  if (!merit)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_vec_free(&merit->sparse_cache));

  SLEQP_CALL(sleqp_vec_free(&merit->multipliers));

  SLEQP_CALL(sleqp_vec_free(&merit->combined_cons_val));

  SLEQP_CALL(sleqp_vec_free(&merit->jac_dot_sparse));

  sleqp_free(&merit->dense_cache);

  SLEQP_CALL(sleqp_params_release(&merit->params));

  SLEQP_CALL(sleqp_problem_release(&merit->problem));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_merit_capture(SleqpMerit* merit)
{
  ++merit->refcount;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_merit_release(SleqpMerit** star)
{
  SleqpMerit* merit = *star;

  if (!merit)
  {
    return SLEQP_OKAY;
  }

  if (--merit->refcount == 0)
  {
    SLEQP_CALL(merit_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
