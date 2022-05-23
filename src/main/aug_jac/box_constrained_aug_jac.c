#include "box_constrained_aug_jac.h"

#include <assert.h>

#include "mem.h"
#include "working_set.h"

typedef struct
{
  SleqpIterate* iterate;
} AugJacData;

static SLEQP_RETCODE
box_constrained_set_iterate(SleqpIterate* iterate, void* data)
{
  AugJacData* jacobian = (AugJacData*)data;

  SLEQP_CALL(sleqp_iterate_release(&jacobian->iterate));

  SLEQP_CALL(sleqp_iterate_capture(iterate));

  jacobian->iterate = iterate;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
box_constrained_min_norm_solution(const SleqpVec* rhs,
                                  SleqpVec* sol,
                                  void* data)
{
  AugJacData* jacobian         = (AugJacData*)data;
  SleqpIterate* iterate        = jacobian->iterate;
  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);

  SLEQP_CALL(sleqp_vec_clear(sol));
  SLEQP_CALL(sleqp_vec_reserve(sol, rhs->nnz));

  for (int k = 0; k < rhs->nnz; ++k)
  {
    const int i        = rhs->indices[k];
    const double value = rhs->data[k];

    const int j = sleqp_working_set_content(working_set, i);

    SLEQP_CALL(sleqp_vec_push(sol, j, value));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
box_constrained_projection(const SleqpVec* rhs,
                           SleqpVec* primal_sol,
                           SleqpVec* dual_sol,
                           void* data)
{
  AugJacData* jacobian         = (AugJacData*)data;
  SleqpIterate* iterate        = jacobian->iterate;
  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);

  if (primal_sol)
  {
    SLEQP_CALL(sleqp_vec_clear(primal_sol));
    SLEQP_CALL(sleqp_vec_reserve(primal_sol, rhs->nnz));
  }

  if (dual_sol)
  {
    SLEQP_CALL(sleqp_vec_clear(dual_sol));
    SLEQP_CALL(sleqp_vec_reserve(dual_sol, rhs->nnz));
  }

  for (int k = 0; k < rhs->nnz; ++k)
  {
    const int j        = rhs->indices[k];
    const double value = rhs->data[k];
    const int index    = sleqp_working_set_var_index(working_set, j);

    if (index == SLEQP_NONE)
    {
      if (primal_sol)
      {
        SLEQP_CALL(sleqp_vec_push(primal_sol, j, value));
      }
    }
    else
    {
      if (dual_sol)
      {
        SLEQP_CALL(sleqp_vec_push(dual_sol, index, value));
      }
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
box_constrained_condition(bool* exact, double* condition, void* data)
{
  *exact     = true;
  *condition = 1.;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
box_constrained_free(void* data)
{
  AugJacData* jacobian = (AugJacData*)data;

  SLEQP_CALL(sleqp_iterate_release(&jacobian->iterate));

  sleqp_free(&jacobian);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
box_constrained_data_create(AugJacData** star)
{
  SLEQP_CALL(sleqp_malloc(star));

  AugJacData* jacobian = *star;

  *jacobian = (AugJacData){0};

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_box_constrained_aug_jac_create(SleqpAugJac** star, SleqpProblem* problem)
{
  AugJacData* aug_jac_data;

  SLEQP_CALL(box_constrained_data_create(&aug_jac_data));

  SleqpAugJacCallbacks callbacks
    = {.set_iterate       = box_constrained_set_iterate,
       .min_norm_solution = box_constrained_min_norm_solution,
       .projection        = box_constrained_projection,
       .condition         = box_constrained_condition,
       .free              = box_constrained_free};

  SLEQP_CALL(
    sleqp_aug_jac_create(star, problem, &callbacks, (void*)aug_jac_data));

  return SLEQP_OKAY;
}
