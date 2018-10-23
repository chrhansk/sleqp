#include "sleqp_soc.h"

#include "sleqp_mem.h"

struct SleqpSOCData
{
  SleqpSparseVec* rhs;
};

SLEQP_RETCODE sleqp_soc_data_create(SleqpSOCData** star,
                                    SleqpProblem* problem)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpSOCData* soc_data = *star;

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

  SLEQP_CALL(sleqp_sparse_vector_clear(rhs));

  int active_set_size = sleqp_aug_jacobian_active_set_size(aug_jacobian);

  SLEQP_CALL(sleqp_sparse_vector_resize(rhs, active_set_size));

  int max_num_nnz = trial_iterate->x->nnz + trial_iterate->cons_val->nnz;

  SLEQP_CALL(sleqp_sparse_vector_reserve(rhs, max_num_nnz));

  {
    SleqpSparseVec* x = trial_iterate->x;

    for(int k_x = 0; k_x < x->nnz; ++k_x)
    {
      int index = x->indices[k_x];

      int variable_index = sleqp_aug_jacobian_variable_index(aug_jacobian, index);

      if(variable_index == -1)
      {
        continue;
      }

      SLEQP_CALL(sleqp_sparse_vector_push(rhs, variable_index, -(x->data[k_x])));

    }
  }


  {
    SleqpSparseVec* cons_val = trial_iterate->cons_val;

    for(int k_c = 0; k_c < cons_val->nnz; ++k_c)
    {
      int index = cons_val->indices[k_c];

      int constraint_index = sleqp_aug_jacobian_constraint_index(aug_jacobian, index);

      if(constraint_index == -1)
      {
        continue;
      }

      SLEQP_CALL(sleqp_sparse_vector_push(rhs, constraint_index, -(cons_val->data[k_c])));
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
