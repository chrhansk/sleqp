#include "dual_estimation_lsq.h"

#include <assert.h>

#include "aug_jac/aug_jac.h"
#include "mem.h"

typedef struct
{
  SleqpProblem* problem;

  SleqpAugJac* aug_jac;

  SleqpVec* solution;
  SleqpVec* neg_grad;

} EstimationData;

static SLEQP_RETCODE
estimate_duals_internal(const SleqpIterate* iterate,
                        SleqpVec* cons_dual,
                        SleqpVec* vars_dual,
                        void* data,
                        int* num_clipped_vars,
                        int* num_clipped_cons)
{
  EstimationData* estimation_data = (EstimationData*)data;
  SleqpAugJac* aug_jac            = estimation_data->aug_jac;

  SleqpProblem* problem        = estimation_data->problem;
  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);

  const int num_variables = sleqp_problem_num_vars(problem);

  SleqpVec* grad = sleqp_iterate_obj_grad(iterate);

  SleqpVec* dual_sol = estimation_data->solution;
  SleqpVec* neg_grad = estimation_data->neg_grad;

  SLEQP_CALL(sleqp_vec_resize(dual_sol, sleqp_working_set_size(working_set)));

  SLEQP_CALL(sleqp_vec_copy(grad, neg_grad));
  SLEQP_CALL(sleqp_vec_scale(neg_grad, -1.));

  SLEQP_CALL(sleqp_aug_jac_solve_lsq(aug_jac, neg_grad, dual_sol));

  *num_clipped_vars = 0;
  *num_clipped_cons = 0;

  {
    SLEQP_CALL(sleqp_vec_reserve(cons_dual, dual_sol->nnz));
    SLEQP_CALL(sleqp_vec_reserve(vars_dual, dual_sol->nnz));

    SLEQP_CALL(sleqp_vec_clear(cons_dual));
    SLEQP_CALL(sleqp_vec_clear(vars_dual));

    for (int k = 0; k < dual_sol->nnz; ++k)
    {
      int sol_index     = dual_sol->indices[k];
      double dual_value = dual_sol->data[k];

      int index = sleqp_working_set_content(working_set, sol_index);

      if (index < num_variables)
      {
        SLEQP_ACTIVE_STATE var_state
          = sleqp_working_set_var_state(working_set, index);

        assert(var_state != SLEQP_INACTIVE);

        if (var_state == SLEQP_ACTIVE_UPPER)
        {
          if (dual_value > 0.)
          {
            SLEQP_CALL(sleqp_vec_push(vars_dual, index, dual_value));
          }
          else if (dual_value < 0.)
          {
            ++(*num_clipped_vars);
          }
        }
        else if (var_state == SLEQP_ACTIVE_LOWER)
        {
          if (dual_value < 0.)
          {
            SLEQP_CALL(sleqp_vec_push(vars_dual, index, dual_value));
          }
          else if (dual_value > 0.)
          {
            ++(*num_clipped_vars);
          }
        }
        else
        {
          assert(var_state == SLEQP_ACTIVE_BOTH);

          SLEQP_CALL(sleqp_vec_push(vars_dual, index, dual_value));
        }
      }
      else
      {
        index -= num_variables;

        SLEQP_ACTIVE_STATE cons_state
          = sleqp_working_set_cons_state(working_set, index);

        assert(cons_state != SLEQP_INACTIVE);

        if (cons_state == SLEQP_ACTIVE_UPPER)
        {
          if (dual_value > 0.)
          {
            SLEQP_CALL(sleqp_vec_push(cons_dual, index, dual_value));
          }
          else if (dual_value < 0.)
          {
            ++(*num_clipped_cons);
          }
        }
        else if (cons_state == SLEQP_ACTIVE_LOWER)
        {
          if (dual_value < 0.)
          {
            SLEQP_CALL(sleqp_vec_push(cons_dual, index, dual_value));
          }
          else if (dual_value > 0.)
          {
            ++(*num_clipped_cons);
          }
        }
        else
        {
          assert(cons_state == SLEQP_ACTIVE_BOTH);

          SLEQP_CALL(sleqp_vec_push(cons_dual, index, dual_value));
        }
      }
    }
  }

  sleqp_log_debug("Dual estimation clipped %d variable and %d constraint duals",
                  *num_clipped_vars,
                  *num_clipped_cons);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
estimate_duals_lsq(const SleqpIterate* iterate,
                   SleqpVec* cons_dual,
                   SleqpVec* vars_dual,
                   void* data)
{
  int num_clipped_vars, num_clipped_cons;

  return estimate_duals_internal(iterate,
                                 cons_dual,
                                 vars_dual,
                                 data,
                                 &num_clipped_vars,
                                 &num_clipped_cons);
}

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_estimate_duals_lsq(SleqpDualEstimation* estimation,
                         const SleqpIterate* iterate,
                         SleqpVec* cons_dual,
                         SleqpVec* vars_dual,
                         int* num_clipped_vars,
                         int* num_clipped_cons)
{
  void* data = sleqp_dual_estimation_data(estimation);

  return estimate_duals_internal(iterate,
                                 cons_dual,
                                 vars_dual,
                                 data,
                                 num_clipped_vars,
                                 num_clipped_cons);
}

static SLEQP_RETCODE
estimation_free(void* data)
{
  EstimationData* estimation_data = (EstimationData*)data;

  SLEQP_CALL(sleqp_vec_free(&estimation_data->solution));
  SLEQP_CALL(sleqp_vec_free(&estimation_data->neg_grad));

  SLEQP_CALL(sleqp_aug_jac_release(&estimation_data->aug_jac));

  SLEQP_CALL(sleqp_problem_release(&estimation_data->problem));

  sleqp_free(&data);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
estimation_data_create(EstimationData** star,
                       SleqpProblem* problem,
                       SleqpAugJac* aug_jac)
{
  SLEQP_CALL(sleqp_malloc(star));

  const int num_vars = sleqp_problem_num_vars(problem);

  EstimationData* estimation_data = *star;

  *estimation_data = (EstimationData){0};

  SLEQP_CALL(sleqp_vec_create_empty(&estimation_data->solution, 0));

  SLEQP_CALL(sleqp_vec_create_empty(&estimation_data->neg_grad, num_vars));

  SLEQP_CALL(sleqp_problem_capture(problem));
  estimation_data->problem = problem;

  SLEQP_CALL(sleqp_aug_jac_capture(aug_jac));
  estimation_data->aug_jac = aug_jac;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_dual_estimation_lsq_create(SleqpDualEstimation** star,
                                 SleqpProblem* problem,
                                 SleqpAugJac* aug_jacobian)
{
  SleqpDualEstimationCallbacks callbacks = {
    .estimate_duals  = estimate_duals_lsq,
    .estimation_free = estimation_free,
  };

  EstimationData* estimation_data;

  SLEQP_CALL(estimation_data_create(&estimation_data, problem, aug_jacobian));

  SLEQP_CALL(
    sleqp_dual_estimation_create(star, &callbacks, (void*)estimation_data));

  return SLEQP_OKAY;
}
