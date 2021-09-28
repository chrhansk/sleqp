#include "unconstrained_aug_jac.h"

#include <assert.h>

static SLEQP_RETCODE
aug_jac_set_iterate(SleqpIterate* iterate,
                    void* aug_jac)
{

#ifndef DEBUG

  SleqpWorkingSet* working_set = sleqp_iterate_get_working_set(iterate);

  assert(sleqp_working_set_size(working_set) == 0);

#endif

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
aug_jac_min_norm_solution(SleqpSparseVec* rhs,
                          SleqpSparseVec* sol,
                          void* aug_jac)
{
  assert(rhs->dim == 0);

  SLEQP_CALL(sleqp_sparse_vector_clear(sol));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
aug_jac_projection(SleqpSparseVec* rhs,
                   SleqpSparseVec* primal_sol,
                   SleqpSparseVec* dual_sol,
                   void* aug_jac)
{

#ifndef DEBUG

  if(dual_sol)
  {
    assert(dual_sol->dim == 0);
  }

#endif

  SLEQP_CALL(sleqp_sparse_vector_copy(rhs, primal_sol));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
aug_jac_condition(bool* exact,
                  double* condition,
                  void *aug_jac)
{
  *exact = true;
  *condition = 1.;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_unconstrained_aug_jac_create(SleqpAugJac** star,
                                                 SleqpProblem* problem)
{
  SleqpAugJacCallbacks callbacks = {
    .set_iterate       = aug_jac_set_iterate,
    .min_norm_solution = aug_jac_min_norm_solution,
    .projection        = aug_jac_projection,
    .condition         = aug_jac_condition,
    .free              = NULL
  };

  SLEQP_CALL(sleqp_aug_jac_create(star,
                                  problem,
                                  &callbacks,
                                  NULL));

  return SLEQP_OKAY;
}
