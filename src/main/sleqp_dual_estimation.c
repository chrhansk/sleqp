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
  SleqpWorkingSet* working_set = sleqp_iterate_get_working_set(iterate);

  SleqpSparseVec* grad = sleqp_iterate_get_func_grad(iterate);

  SleqpSparseVec* dual_sol = estimation_data->solution;
  SleqpSparseVec* neg_grad = estimation_data->neg_grad;

  SLEQP_CALL(sleqp_sparse_vector_resize(dual_sol,
                                        sleqp_working_set_size(working_set)));

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

  int num_clipped_vars = 0, num_clipped_cons = 0;

  {
    SleqpSparseVec* cons_dual = sleqp_iterate_get_cons_dual(iterate);
    SleqpSparseVec* vars_dual = sleqp_iterate_get_vars_dual(iterate);

    SLEQP_CALL(sleqp_sparse_vector_reserve(cons_dual, dual_sol->nnz));
    SLEQP_CALL(sleqp_sparse_vector_reserve(vars_dual, dual_sol->nnz));

    vars_dual->nnz = 0;
    cons_dual->nnz = 0;

    for(int k = 0; k < dual_sol->nnz; ++k)
    {
      int sol_index = dual_sol->indices[k];
      double dual_value = dual_sol->data[k];

      int index = sleqp_working_set_get_content(working_set, sol_index);

      if(index < problem->num_variables)
      {
        SLEQP_ACTIVE_STATE var_state = sleqp_working_set_get_variable_state(working_set, index);

        assert(var_state != SLEQP_INACTIVE);

        if(var_state == SLEQP_ACTIVE_UPPER)
        {
          if(dual_value > 0.)
          {
            SLEQP_CALL(sleqp_sparse_vector_push(vars_dual, index, dual_value));
          }
          else if(dual_value < 0.)
          {
            ++num_clipped_vars;
          }
        }
        else if(var_state == SLEQP_ACTIVE_LOWER)
        {
          if(dual_value < 0.)
          {
            SLEQP_CALL(sleqp_sparse_vector_push(vars_dual, index, dual_value));
          }
          else if(dual_value > 0.)
          {
            ++num_clipped_vars;
          }
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

        SLEQP_ACTIVE_STATE cons_state = sleqp_working_set_get_constraint_state(working_set, index);

        assert(cons_state != SLEQP_INACTIVE);

        if(cons_state == SLEQP_ACTIVE_UPPER)
        {
          if(dual_value > 0.)
          {
            SLEQP_CALL(sleqp_sparse_vector_push(cons_dual, index, dual_value));
          }
          else if(dual_value < 0.)
          {
            ++num_clipped_cons;
          }
        }
        else if(cons_state == SLEQP_ACTIVE_LOWER)
        {
          if(dual_value < 0.)
          {
            SLEQP_CALL(sleqp_sparse_vector_push(cons_dual, index, dual_value));
          }
          else if(dual_value > 0.)
          {
            ++num_clipped_cons;
          }
        }
        else
        {
          assert(cons_state == SLEQP_ACTIVE_BOTH);

          SLEQP_CALL(sleqp_sparse_vector_push(cons_dual, index, dual_value));
        }
      }
    }
  }

  sleqp_log_debug("Dual estimation clipped %d variable and %d constraint duals",
                  num_clipped_vars,
                  num_clipped_cons);

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
