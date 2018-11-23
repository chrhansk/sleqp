#include "sleqp_soc.h"

#include "sleqp_cmp.h"
#include "sleqp_mem.h"

struct SleqpSOCData
{
  SleqpProblem* problem;
  SleqpParams* params;

  SleqpSparseVec* upper_diff;
  SleqpSparseVec* lower_diff;

  SleqpSparseVec* rhs;

};

SLEQP_RETCODE sleqp_soc_data_create(SleqpSOCData** star,
                                    SleqpProblem* problem,
                                    SleqpParams* params)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpSOCData* soc_data = *star;

  soc_data->problem = problem;
  soc_data->params = params;

  SLEQP_CALL(sleqp_sparse_vector_create(&soc_data->upper_diff, 0, 0));
  SLEQP_CALL(sleqp_sparse_vector_create(&soc_data->lower_diff, 0, 0));

  SLEQP_CALL(sleqp_sparse_vector_create(&soc_data->rhs, 0, 0));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE add_variable_entries(SleqpSOCData* soc_data,
                                          SleqpIterate* iterate,
                                          SleqpIterate* trial_iterate)
{
  SleqpProblem* problem = soc_data->problem;
  const double eps = sleqp_params_get_eps(soc_data->params);

  SleqpSparseVec* lower_diff = soc_data->lower_diff;
  SleqpSparseVec* upper_diff = soc_data->upper_diff;

  SleqpSparseVec* rhs = soc_data->rhs;
  SleqpActiveSet* active_set = iterate->active_set;

  SLEQP_CALL(sleqp_sparse_vector_clear(soc_data->lower_diff));
  SLEQP_CALL(sleqp_sparse_vector_resize(soc_data->lower_diff, problem->num_variables));

  SLEQP_CALL(sleqp_sparse_vector_clear(soc_data->upper_diff));
  SLEQP_CALL(sleqp_sparse_vector_resize(soc_data->upper_diff, problem->num_variables));

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(trial_iterate->x,
                                            problem->var_lb,
                                            -1.,
                                            1.,
                                            eps,
                                            lower_diff));

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(trial_iterate->x,
                                            problem->var_ub,
                                            -1.,
                                            1.,
                                            eps,
                                            upper_diff));

  int k_ld = 0., k_ud = 0.;

  for(int i = 0; i < problem->num_variables; ++i)
  {
    while(k_ld < lower_diff->nnz && lower_diff->indices[k_ld] < i)
    {
      ++k_ld;
    }

    while(k_ud < upper_diff->nnz && upper_diff->indices[k_ud] < i)
    {
      ++k_ud;
    }

    bool valid_ud = (k_ud < upper_diff->nnz && upper_diff->indices[k_ud] == i);
    bool valid_ld = (k_ld < lower_diff->nnz && lower_diff->indices[k_ld] == i);

    const double udval = valid_ud ? upper_diff->data[k_ud] : 0.;
    const double ldval = valid_ld ? lower_diff->data[k_ld] : 0.;

    const SLEQP_ACTIVE_STATE variable_state = sleqp_active_set_get_variable_state(active_set, i);

    if(variable_state == SLEQP_INACTIVE)
    {
      continue;
    }

    const int variable_index = sleqp_active_set_get_variable_index(active_set, i);

    if(variable_state == SLEQP_ACTIVE_UPPER && !sleqp_zero(udval, eps))
    {
      SLEQP_CALL(sleqp_sparse_vector_push(rhs, variable_index, udval));
    }
    else if(variable_state == SLEQP_ACTIVE_LOWER && !sleqp_zero(ldval, eps))
    {
      SLEQP_CALL(sleqp_sparse_vector_push(rhs, variable_index, -ldval));
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE add_constraint_entries(SleqpSOCData* soc_data,
                                            SleqpIterate* iterate,
                                            SleqpIterate* trial_iterate)
{
  SleqpProblem* problem = soc_data->problem;
  const double eps = sleqp_params_get_eps(soc_data->params);

  SleqpSparseVec* lower_diff = soc_data->lower_diff;
  SleqpSparseVec* upper_diff = soc_data->upper_diff;

  SleqpSparseVec* rhs = soc_data->rhs;
  SleqpActiveSet* active_set = iterate->active_set;

  SLEQP_CALL(sleqp_sparse_vector_clear(soc_data->lower_diff));
  SLEQP_CALL(sleqp_sparse_vector_resize(soc_data->lower_diff, problem->num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_clear(soc_data->upper_diff));
  SLEQP_CALL(sleqp_sparse_vector_resize(soc_data->upper_diff, problem->num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(trial_iterate->cons_val,
                                            problem->cons_lb,
                                            -1.,
                                            1.,
                                            eps,
                                            lower_diff));

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(trial_iterate->cons_val,
                                            problem->cons_ub,
                                            -1.,
                                            1.,
                                            eps,
                                            upper_diff));


  int k_ld = 0., k_ud = 0.;

  for(int i = 0; i < problem->num_constraints; ++i)
  {
    while(k_ld < lower_diff->nnz && lower_diff->indices[k_ld] < i)
    {
      ++k_ld;
    }

    while(k_ud < upper_diff->nnz && upper_diff->indices[k_ud] < i)
    {
      ++k_ud;
    }

    bool valid_ud = (k_ud < upper_diff->nnz && upper_diff->indices[k_ud] == i);
    bool valid_ld = (k_ld < lower_diff->nnz && lower_diff->indices[k_ld] == i);

    const double udval = valid_ud ? upper_diff->data[k_ud] : 0.;
    const double ldval = valid_ld ? lower_diff->data[k_ld] : 0.;

    const SLEQP_ACTIVE_STATE constraint_state = sleqp_active_set_get_constraint_state(active_set, i);

    if(constraint_state == SLEQP_INACTIVE)
    {
      continue;
    }

    const int constraint_index = sleqp_active_set_get_constraint_index(active_set, i);

    if(constraint_state == SLEQP_ACTIVE_UPPER && !sleqp_zero(udval, eps))
    {
      SLEQP_CALL(sleqp_sparse_vector_push(rhs, constraint_index, udval));
    }
    else if(constraint_state == SLEQP_ACTIVE_LOWER && !sleqp_zero(ldval, eps))
    {
      SLEQP_CALL(sleqp_sparse_vector_push(rhs, constraint_index, -ldval));
    }
  }

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

  SLEQP_CALL(add_variable_entries(soc_data, iterate, trial_iterate));
  SLEQP_CALL(add_constraint_entries(soc_data, iterate, trial_iterate));

  SLEQP_CALL(sleqp_aug_jacobian_min_norm_solution(aug_jacobian, rhs, soc_direction));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_soc_data_free(SleqpSOCData** star)
{
  SleqpSOCData* soc_data = *star;

  SLEQP_CALL(sleqp_sparse_vector_free(&soc_data->lower_diff));
  SLEQP_CALL(sleqp_sparse_vector_free(&soc_data->upper_diff));

  SLEQP_CALL(sleqp_sparse_vector_free(&soc_data->rhs));

  sleqp_free(star);

  return SLEQP_OKAY;
}
