#include "box_constrained_aug_jac.h"

#include <assert.h>

#include "mem.h"
#include "working_set.h"

typedef struct
{
  SleqpIterate* iterate;
} AugJacData;

static SLEQP_RETCODE
aug_jac_set_iterate(SleqpIterate* iterate,
                    void* data)
{
  AugJacData* jacobian = (AugJacData*) data;

  SLEQP_CALL(sleqp_iterate_release(&jacobian->iterate));

  SLEQP_CALL(sleqp_iterate_capture(iterate));

  jacobian->iterate = iterate;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
aug_jac_min_norm_solution(SleqpSparseVec* rhs,
                          SleqpSparseVec* sol,
                          void* data)
{
  AugJacData* jacobian = (AugJacData*) data;
  SleqpIterate* iterate = jacobian->iterate;
  SleqpWorkingSet* working_set = sleqp_iterate_get_working_set(iterate);

  SLEQP_CALL(sleqp_sparse_vector_clear(sol));

  for(int k = 0; k < rhs->nnz ; ++k)
  {
    const int j = rhs->indices[k];
    const double value = rhs->data[k];
    const int index = sleqp_working_set_get_variable_index(working_set, j);

    if(index == SLEQP_NONE)
    {
      SLEQP_CALL(sleqp_sparse_vector_push(sol, j, value));
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
aug_jac_projection(SleqpSparseVec* rhs,
                   SleqpSparseVec* primal_sol,
                   SleqpSparseVec* dual_sol,
                   void* data)
{
  AugJacData* jacobian = (AugJacData*) data;
  SleqpIterate* iterate = jacobian->iterate;
  SleqpWorkingSet* working_set = sleqp_iterate_get_working_set(iterate);

  if(primal_sol)
  {
    SLEQP_CALL(sleqp_sparse_vector_clear(primal_sol));
    SLEQP_CALL(sleqp_sparse_vector_reserve(primal_sol, rhs->nnz));
  }

  if(dual_sol)
  {
    SLEQP_CALL(sleqp_sparse_vector_clear(dual_sol));
    SLEQP_CALL(sleqp_sparse_vector_reserve(dual_sol, rhs->nnz));
  }

  for(int k = 0; k < rhs->nnz ; ++k)
  {
    const int j = rhs->indices[k];
    const double value = rhs->data[k];
    const int index = sleqp_working_set_get_variable_index(working_set, j);

    if(index == SLEQP_NONE)
    {
      if(primal_sol)
      {
        SLEQP_CALL(sleqp_sparse_vector_push(primal_sol, j, value));
      }
    }
    else
    {
      if(dual_sol)
      {
        SLEQP_CALL(sleqp_sparse_vector_push(dual_sol, index, value));
      }
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
aug_jac_condition(bool* exact,
                  double* condition,
                  void* data)
{
  *exact = true;
  *condition = 1.;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
aug_jac_free(void* data)
{
  AugJacData* jacobian = (AugJacData*) data;

  SLEQP_CALL(sleqp_iterate_release(&jacobian->iterate));

  sleqp_free(&jacobian);

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE aug_jac_data_create(AugJacData** star)
{
  SLEQP_CALL(sleqp_malloc(star));

  AugJacData* jacobian = *star;

  *jacobian = (AugJacData){0};

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_box_constrained_aug_jac_create(SleqpAugJac** star,
                                                   SleqpProblem* problem)
{
  AugJacData* aug_jac_data;

  SLEQP_CALL(aug_jac_data_create(&aug_jac_data));

  SleqpAugJacCallbacks callbacks = {
    .set_iterate       = aug_jac_set_iterate,
    .min_norm_solution = aug_jac_min_norm_solution,
    .projection        = aug_jac_projection,
    .condition         = aug_jac_condition,
    .free              = aug_jac_free
  };

  SLEQP_CALL(sleqp_aug_jac_create(star,
                                  problem,
                                  &callbacks,
                                  (void*) aug_jac_data));

  return SLEQP_OKAY;
}
