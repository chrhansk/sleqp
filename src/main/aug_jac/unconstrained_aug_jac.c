#include "unconstrained_aug_jac.h"

#include "fail.h"

static SLEQP_RETCODE
aug_jac_set_iterate(SleqpIterate* iterate, void* aug_jac)
{

#ifndef DEBUG

  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);

  assert(sleqp_working_set_size(working_set) == 0);

#endif

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
aug_jac_solve_min_norm(const SleqpVec* rhs, SleqpVec* sol, void* aug_jac)
{
  assert(rhs->dim == 0);
  SLEQP_CALL(sleqp_vec_clear(sol));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
aug_jac_project_nullspace(const SleqpVec* rhs, SleqpVec* sol, void* aug_jac)
{
  SLEQP_CALL(sleqp_vec_copy(rhs, sol));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
aug_jac_solve_lsq(const SleqpVec* rhs, SleqpVec* sol, void* aug_jac)
{
  assert(sol->dim == 0);
  return SLEQP_OKAY;
}

static SLEQP_RETCODE
aug_jac_condition(bool* exact, double* condition, void* aug_jac)
{
  *exact     = true;
  *condition = 1.;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_unconstrained_aug_jac_create(SleqpAugJac** star, SleqpProblem* problem)
{
  SleqpAugJacCallbacks callbacks
    = {.set_iterate       = aug_jac_set_iterate,
       .solve_min_norm    = aug_jac_solve_min_norm,
       .solve_lsq         = aug_jac_solve_lsq,
       .project_nullspace = aug_jac_project_nullspace,
       .condition         = aug_jac_condition,
       .free              = NULL};

  SLEQP_CALL(sleqp_aug_jac_create(star, problem, &callbacks, NULL));

  return SLEQP_OKAY;
}
