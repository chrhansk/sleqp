#include "aug_jac.h"

#include "fail.h"
#include "mem.h"

struct SleqpAugJac
{
  int refcount;

  SleqpProblem* problem;
  int working_set_size;
  SleqpTimer* creation_timer;
  SleqpTimer* solution_timer;

  SleqpAugJacCallbacks callbacks;
  void* data;
};

SLEQP_RETCODE
sleqp_aug_jac_create(SleqpAugJac** star,
                     SleqpProblem* problem,
                     SleqpAugJacCallbacks* callbacks,
                     void* aug_jac_data)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpAugJac* aug_jac = *star;

  *aug_jac = (SleqpAugJac){0};

  aug_jac->refcount = 1;

  SLEQP_CALL(sleqp_problem_capture(problem));
  aug_jac->problem = problem;

  SLEQP_CALL(sleqp_timer_create(&aug_jac->creation_timer));
  SLEQP_CALL(sleqp_timer_create(&aug_jac->solution_timer));

  aug_jac->callbacks = *callbacks;
  aug_jac->data      = aug_jac_data;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_aug_jac_set_iterate(SleqpAugJac* aug_jac, SleqpIterate* iterate)
{
  SLEQP_CALL(sleqp_timer_start(aug_jac->creation_timer));

  SLEQP_CALL(aug_jac->callbacks.set_iterate(iterate, aug_jac->data));

  SLEQP_CALL(sleqp_timer_stop(aug_jac->creation_timer));

  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);

  aug_jac->working_set_size = sleqp_working_set_size(working_set);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_aug_jac_solve_min_norm(SleqpAugJac* aug_jac,
                             const SleqpVec* rhs,
                             SleqpVec* sol)
{
  assert(sol->dim == sleqp_problem_num_vars(aug_jac->problem));
  assert(rhs->dim == aug_jac->working_set_size);

  SLEQP_CALL(sleqp_timer_start(aug_jac->solution_timer));

  SLEQP_CALL(aug_jac->callbacks.solve_min_norm(rhs, sol, aug_jac->data));

  SLEQP_CALL(sleqp_timer_stop(aug_jac->solution_timer));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_aug_jac_solve_lsq(SleqpAugJac* aug_jac,
                        const SleqpVec* rhs,
                        SleqpVec* sol)
{
  assert(rhs->dim == sleqp_problem_num_vars(aug_jac->problem));
  assert(sol->dim == aug_jac->working_set_size);

  return aug_jac->callbacks.solve_lsq(rhs, sol, aug_jac->data);
}

SLEQP_RETCODE
sleqp_aug_jac_project_nullspace(SleqpAugJac* aug_jac,
                                const SleqpVec* rhs,
                                SleqpVec* sol)
{
  assert(rhs->dim == sleqp_problem_num_vars(aug_jac->problem));
  assert(sol->dim == sleqp_problem_num_vars(aug_jac->problem));

  return aug_jac->callbacks.project_nullspace(rhs, sol, aug_jac->data);
}

SLEQP_RETCODE
sleqp_aug_jac_condition(SleqpAugJac* aug_jac, bool* exact, double* condition)
{
  SLEQP_CALL(aug_jac->callbacks.condition(exact, condition, aug_jac->data));

  return SLEQP_OKAY;
}

SleqpTimer*
sleqp_aug_jac_creation_timer(SleqpAugJac* aug_jac)
{
  return aug_jac->creation_timer;
}

SleqpTimer*
sleqp_aug_jac_solution_timer(SleqpAugJac* aug_jac)
{
  return aug_jac->solution_timer;
}

SLEQP_RETCODE
sleqp_aug_jac_capture(SleqpAugJac* aug_jac)
{
  ++aug_jac->refcount;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
aug_jac_free(SleqpAugJac** star)
{
  SleqpAugJac* aug_jac = *star;

  if (aug_jac->callbacks.free)
  {
    SLEQP_CALL(aug_jac->callbacks.free(aug_jac->data));
  }

  SLEQP_CALL(sleqp_problem_release(&aug_jac->problem));

  SLEQP_CALL(sleqp_timer_free(&aug_jac->solution_timer));
  SLEQP_CALL(sleqp_timer_free(&aug_jac->creation_timer));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_aug_jac_release(SleqpAugJac** star)
{
  SleqpAugJac* aug_jac = *star;

  if (!aug_jac)
  {
    return SLEQP_OKAY;
  }

  if (--aug_jac->refcount == 0)
  {
    SLEQP_CALL(aug_jac_free(star));
  }

  return SLEQP_OKAY;
}
