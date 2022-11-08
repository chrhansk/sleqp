#include "direction.h"

#include "cmp.h"
#include "mem.h"
#include "problem.h"
#include "pub_iterate.h"
#include "pub_problem.h"
#include "pub_types.h"
#include "sparse/mat.h"
#include "sparse/pub_vec.h"

struct SleqpDirection
{
  int refcount;
  SleqpVec* primal;
  double obj_grad;
  SleqpVec* cons_jac;
  SleqpVec* hess;
};

SLEQP_RETCODE
sleqp_direction_create(SleqpDirection** star,
                       SleqpProblem* problem,
                       SleqpParams* params)
{
  SLEQP_CALL(sleqp_malloc(star));

  const int num_variables   = sleqp_problem_num_vars(problem);
  const int num_constraints = sleqp_problem_num_cons(problem);

  SleqpDirection* direction = *star;

  *direction = (SleqpDirection){0};

  direction->refcount = 1;

  SLEQP_CALL(sleqp_vec_create_empty(&direction->primal, num_variables));
  SLEQP_CALL(sleqp_vec_create_empty(&direction->cons_jac, num_constraints));
  SLEQP_CALL(sleqp_vec_create_empty(&direction->hess, num_variables));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_direction_reset(SleqpDirection* direction,
                      SleqpProblem* problem,
                      const SleqpIterate* iterate,
                      const SleqpVec* multipliers,
                      double* cache,
                      double zero_eps)
{
  SleqpVec* direction_primal   = sleqp_direction_primal(direction);
  SleqpVec* direction_cons_jac = sleqp_direction_cons_jac(direction);
  SleqpVec* direction_hess     = sleqp_direction_hess(direction);

  const int num_cons = sleqp_problem_num_cons(problem);

  SLEQP_CALL(sleqp_vec_dot(direction_primal,
                           sleqp_iterate_obj_grad(iterate),
                           sleqp_direction_obj_grad(direction)));

  if (num_cons > 0)
  {
    SleqpMat* cons_jac = sleqp_iterate_cons_jac(iterate);

    SLEQP_CALL(sleqp_mat_mult_vec(cons_jac, direction_primal, cache));

    SLEQP_CALL(
      sleqp_vec_set_from_raw(direction_cons_jac, cache, num_cons, zero_eps));
  }
  else
  {
    SLEQP_CALL(sleqp_vec_clear(direction_cons_jac));
  }

  const double one = 1.;
  SLEQP_CALL(sleqp_problem_hess_prod(problem,
                                     &one,
                                     direction_primal,
                                     multipliers,
                                     direction_hess));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_direction_check(const SleqpDirection* direction,
                      SleqpProblem* problem,
                      const SleqpIterate* iterate,
                      const SleqpVec* multipliers,
                      double* cache,
                      double zero_eps,
                      bool* valid)
{
  *valid = true;
  double objective_dot;

  const int num_vars = sleqp_problem_num_vars(problem);
  const int num_cons = sleqp_problem_num_cons(problem);

  SLEQP_CALL(sleqp_vec_dot(direction->primal,
                           sleqp_iterate_obj_grad(iterate),
                           &objective_dot));

  if (!sleqp_is_eq(objective_dot, direction->obj_grad, 1e-6))
  {
    *valid = false;
    return SLEQP_OKAY;
  }

  SleqpVec* cons_dir;

  SLEQP_CALL(sleqp_vec_create_full(&cons_dir, num_cons));

  SLEQP_CALL(sleqp_mat_mult_vec(sleqp_iterate_cons_jac(iterate),
                                direction->primal,
                                cache));

  SLEQP_CALL(sleqp_vec_set_from_raw(cons_dir, cache, num_cons, zero_eps));

  if (!sleqp_vec_eq(cons_dir, direction->cons_jac, 1e-6))
  {
    *valid = false;
    return SLEQP_OKAY;
  }

  const double one = 1.;

  SleqpVec* hess_dir;

  SLEQP_CALL(sleqp_vec_create_full(&hess_dir, num_vars));

  SLEQP_CALL(sleqp_problem_hess_prod(problem,
                                     &one,
                                     direction->primal,
                                     multipliers,
                                     hess_dir));

  if (!sleqp_vec_eq(hess_dir, direction->hess, 1e-6))
  {
    *valid = false;
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_vec_free(&hess_dir));

  SLEQP_CALL(sleqp_vec_free(&cons_dir));

  return SLEQP_OKAY;
}

SleqpVec*
sleqp_direction_primal(const SleqpDirection* direction)
{
  return direction->primal;
}

double*
sleqp_direction_obj_grad(const SleqpDirection* direction)
{
  return (double*)&(direction->obj_grad);
}

SleqpVec*
sleqp_direction_cons_jac(const SleqpDirection* direction)
{
  return direction->cons_jac;
}

SleqpVec*
sleqp_direction_hess(const SleqpDirection* direction)
{
  return direction->hess;
}

SLEQP_RETCODE
sleqp_direction_scale(SleqpDirection* direction, double factor)
{
  SLEQP_CALL(sleqp_vec_scale(direction->primal, factor));
  SLEQP_CALL(sleqp_vec_scale(direction->cons_jac, factor));
  SLEQP_CALL(sleqp_vec_scale(direction->hess, factor));

  direction->obj_grad *= factor;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_direction_set_zero(SleqpDirection* direction)
{
  SLEQP_CALL(sleqp_vec_clear(direction->primal));
  SLEQP_CALL(sleqp_vec_clear(direction->cons_jac));
  SLEQP_CALL(sleqp_vec_clear(direction->hess));

  direction->obj_grad = 0;

  return SLEQP_OKAY;
}

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_direction_add_scaled(const SleqpDirection* first,
                           const SleqpDirection* second,
                           const double first_factor,
                           const double second_factor,
                           const double eps,
                           SleqpDirection* result)
{
  SLEQP_CALL(sleqp_vec_add_scaled(first->primal,
                                  second->primal,
                                  first_factor,
                                  second_factor,
                                  eps,
                                  result->primal));

  SLEQP_CALL(sleqp_vec_add_scaled(first->cons_jac,
                                  second->cons_jac,
                                  first_factor,
                                  second_factor,
                                  eps,
                                  result->cons_jac));

  SLEQP_CALL(sleqp_vec_add_scaled(first->hess,
                                  second->hess,
                                  first_factor,
                                  second_factor,
                                  eps,
                                  result->hess));

  result->obj_grad
    = (first_factor * first->obj_grad) + (second_factor * second->obj_grad);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_direction_copy(const SleqpDirection* source, SleqpDirection* target)
{
  assert(source != target);

  SLEQP_CALL(sleqp_vec_copy(source->primal, target->primal));
  SLEQP_CALL(sleqp_vec_copy(source->cons_jac, target->cons_jac));
  SLEQP_CALL(sleqp_vec_copy(source->hess, target->hess));

  target->obj_grad = source->obj_grad;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_direction_capture(SleqpDirection* direction)
{
  ++direction->refcount;
  return SLEQP_OKAY;
}

static SLEQP_RETCODE
direction_free(SleqpDirection** star)
{
  SleqpDirection* direction = *star;

  SLEQP_CALL(sleqp_vec_free(&direction->hess));
  SLEQP_CALL(sleqp_vec_free(&direction->cons_jac));
  SLEQP_CALL(sleqp_vec_free(&direction->primal));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_direction_release(SleqpDirection** star)
{
  SleqpDirection* direction = *star;

  if (!direction)
  {
    return SLEQP_OKAY;
  }

  if (--direction->refcount == 0)
  {
    SLEQP_CALL(direction_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
