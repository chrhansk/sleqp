#include "sleqp_dual_estimation.h"

#include <umfpack.h>

#include "sleqp_cmp.h"
#include "sleqp_mem.h"

struct SleqpDualEstimationData
{
  SleqpProblem* problem;

  SleqpSparseVec* solution;
  SleqpSparseVec* neg_grad;
};


SLEQP_RETCODE sleqp_dual_estimation_data_create(SleqpDualEstimationData** star,
                                                SleqpProblem* problem)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpDualEstimationData* data = *star;

  data->problem = problem;

  SLEQP_CALL(sleqp_sparse_vector_create(&data->solution, 0, 0));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->neg_grad, problem->num_variables, 0));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_dual_estimation_compute(SleqpDualEstimationData* estimation_data,
                                            SleqpIterate* iterate,
                                            SleqpSparseVec* residuum,
                                            SleqpAugJacobian* jacobian)
{
  SleqpProblem* problem = estimation_data->problem;
  SleqpActiveSet* active_set = iterate->active_set;

  SleqpSparseVec* grad = iterate->func_grad;

  SleqpSparseVec* dual_sol = estimation_data->solution;
  SleqpSparseVec* neg_grad = estimation_data->neg_grad;

  SLEQP_CALL(sleqp_sparse_vector_resize(dual_sol,
                                        sleqp_active_set_size(iterate->active_set)));

  SLEQP_CALL(sleqp_sparse_vector_clear(neg_grad));

  SLEQP_CALL(sleqp_sparse_vector_reserve(neg_grad, grad->nnz));

  for(int k = 0; k < grad->nnz; ++k)
  {
    SLEQP_CALL(sleqp_sparse_vector_push(neg_grad, grad->indices[k], -1.*grad->data[k]));
  }

  SLEQP_CALL(sleqp_aug_jacobian_projection(jacobian,
                                           neg_grad,
                                           residuum,
                                           dual_sol));

  {
    SleqpSparseVec* cons_dual = iterate->cons_dual;
    SleqpSparseVec* vars_dual = iterate->vars_dual;

    SLEQP_CALL(sleqp_sparse_vector_reserve(cons_dual, dual_sol->nnz));
    SLEQP_CALL(sleqp_sparse_vector_reserve(vars_dual, dual_sol->nnz));

    vars_dual->nnz = 0;
    cons_dual->nnz = 0;

    for(int k = 0; k < dual_sol->nnz; ++k)
    {
      int sol_index = dual_sol->indices[k];
      double dual_value = dual_sol->data[k];

      int index = sleqp_active_set_get_content(active_set, sol_index);

      if(index < problem->num_variables)
      {
        SLEQP_ACTIVE_STATE var_state = sleqp_active_set_get_variable_state(active_set, index);

        assert(var_state != SLEQP_INACTIVE);

        if(var_state == SLEQP_ACTIVE_UPPER)
        {
          SLEQP_CALL(sleqp_sparse_vector_push(vars_dual, index, SLEQP_MAX(dual_value, 0.)));
        }
        else if(var_state == SLEQP_ACTIVE_LOWER)
        {
          SLEQP_CALL(sleqp_sparse_vector_push(vars_dual, index, SLEQP_MIN(dual_value, 0.)));
        }
        else
        {
          assert(var_state == SLEQP_ACTIVE_BOTH);

          SLEQP_CALL(sleqp_sparse_vector_push(vars_dual, index, dual_value));
        }

      }
      else
      {
        index -= problem->num_variables;

        SLEQP_ACTIVE_STATE cons_state = sleqp_active_set_get_constraint_state(active_set, index);

        assert(cons_state != SLEQP_INACTIVE);

        if(cons_state == SLEQP_ACTIVE_UPPER)
        {
          SLEQP_CALL(sleqp_sparse_vector_push(cons_dual, index, SLEQP_MAX(dual_value, 0.)));
        }
        else if(cons_state == SLEQP_ACTIVE_LOWER)
        {
          SLEQP_CALL(sleqp_sparse_vector_push(cons_dual, index, SLEQP_MIN(dual_value, 0.)));
        }
        else
        {
          assert(cons_state == SLEQP_ACTIVE_BOTH);

          SLEQP_CALL(sleqp_sparse_vector_push(cons_dual, index, dual_value));
        }
      }

    }

  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_dual_estimation_data_free(SleqpDualEstimationData** star)
{
  SleqpDualEstimationData* data = *star;

  if(!data)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_sparse_vector_free(&data->neg_grad));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->solution));

  sleqp_free(star);

  return SLEQP_OKAY;
}
