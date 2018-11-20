#include "sleqp_soc.h"

#include "sleqp_mem.h"

struct SleqpSOCData
{
  SleqpProblem* problem;
  SleqpSparseVec* rhs;
};

SLEQP_RETCODE sleqp_soc_data_create(SleqpSOCData** star,
                                    SleqpProblem* problem)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpSOCData* soc_data = *star;

  soc_data->problem = problem;

  SLEQP_CALL(sleqp_sparse_vector_create(&(soc_data->rhs), 0, 0));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_soc_compute(SleqpSOCData* soc_data,
                                SleqpAugJacobian* aug_jacobian,
                                SleqpIterate* iterate,
                                SleqpIterate* trial_iterate,
                                SleqpSparseVec* soc_direction)
{
  SleqpSparseVec* rhs = soc_data->rhs;
  SleqpActiveSet* active_set = iterate->active_set;

  const int active_set_size = sleqp_active_set_size(active_set);

  SLEQP_CALL(sleqp_sparse_vector_clear(rhs));

  SLEQP_CALL(sleqp_sparse_vector_resize(rhs, active_set_size));

  SLEQP_CALL(sleqp_sparse_vector_reserve(rhs, active_set_size));

  int max_num_nnz = trial_iterate->x->nnz + trial_iterate->cons_val->nnz;

  SLEQP_CALL(sleqp_sparse_vector_reserve(rhs, max_num_nnz));

  {
    SleqpSparseVec* x = trial_iterate->x;

    for(int k_x = 0; k_x < x->nnz; ++k_x)
    {
      const int index = x->indices[k_x];

      const SLEQP_ACTIVE_STATE variable_state = sleqp_active_set_get_variable_state(active_set, index);

      if(variable_state == SLEQP_INACTIVE)
      {
        continue;
      }

      const int variable_index = sleqp_active_set_get_variable_index(active_set, index);

      const double val = (variable_state == SLEQP_ACTIVE_UPPER) ? x->data[k_x] : (-1. * x->data[k_x]);

      assert(variable_index != -1);

      SLEQP_CALL(sleqp_sparse_vector_push(rhs, variable_index, -val));
    }
  }


  {
    SleqpSparseVec* cons_val = trial_iterate->cons_val;

    for(int k_c = 0; k_c < cons_val->nnz; ++k_c)
    {
      const int index = cons_val->indices[k_c];

      const SLEQP_ACTIVE_STATE constraint_state = sleqp_active_set_get_constraint_state(active_set, index);

      if(constraint_state == SLEQP_INACTIVE)
      {
        continue;
      }

      int constraint_index = sleqp_active_set_get_constraint_index(active_set, index);

      assert(constraint_index != -1);

      const double val = (constraint_state == SLEQP_ACTIVE_UPPER) ? cons_val->data[k_c] : (-1. * cons_val->data[k_c]);

      SLEQP_CALL(sleqp_sparse_vector_push(rhs, constraint_index, -val));
    }
  }

  SLEQP_CALL(sleqp_aug_jacobian_min_norm_solution(aug_jacobian, rhs, soc_direction));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_soc_data_free(SleqpSOCData** star)
{
  SleqpSOCData* soc_data = *star;

  SLEQP_CALL(sleqp_sparse_vector_free(&(soc_data->rhs)));

  sleqp_free(star);

  return SLEQP_OKAY;
}
