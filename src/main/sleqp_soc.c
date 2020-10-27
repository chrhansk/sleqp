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

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&soc_data->upper_diff, 0));
  SLEQP_CALL(sleqp_sparse_vector_create_empty(&soc_data->lower_diff, 0));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&soc_data->rhs, 0));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE add_variable_entries(SleqpSOCData* soc_data,
                                          SleqpIterate* iterate,
                                          SleqpIterate* trial_iterate)
{
  SleqpProblem* problem = soc_data->problem;

  const double eps = sleqp_params_get_eps(soc_data->params);
  const double zero_eps = sleqp_params_get_zero_eps(soc_data->params);

  SleqpSparseVec* lower_diff = soc_data->lower_diff;
  SleqpSparseVec* upper_diff = soc_data->upper_diff;

  SleqpSparseVec* rhs = soc_data->rhs;
  SleqpWorkingSet* working_set = sleqp_iterate_get_working_set(iterate);

  SLEQP_CALL(sleqp_sparse_vector_clear(soc_data->lower_diff));
  SLEQP_CALL(sleqp_sparse_vector_resize(soc_data->lower_diff, problem->num_variables));

  SLEQP_CALL(sleqp_sparse_vector_clear(soc_data->upper_diff));
  SLEQP_CALL(sleqp_sparse_vector_resize(soc_data->upper_diff, problem->num_variables));

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(sleqp_iterate_get_primal(trial_iterate),
                                            problem->var_lb,
                                            -1.,
                                            1.,
                                            zero_eps,
                                            lower_diff));

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(sleqp_iterate_get_primal(trial_iterate),
                                            problem->var_ub,
                                            -1.,
                                            1.,
                                            zero_eps,
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

    const SLEQP_ACTIVE_STATE variable_state = sleqp_working_set_get_variable_state(working_set, i);

    if(variable_state == SLEQP_INACTIVE)
    {
      continue;
    }

    const int variable_index = sleqp_working_set_get_variable_index(working_set, i);

    if(variable_state == SLEQP_ACTIVE_UPPER && !sleqp_is_zero(udval, eps))
    {
      SLEQP_CALL(sleqp_sparse_vector_push(rhs, variable_index, udval));
    }
    else if(variable_state == SLEQP_ACTIVE_LOWER && !sleqp_is_zero(ldval, eps))
    {
      SLEQP_CALL(sleqp_sparse_vector_push(rhs, variable_index, ldval));
    }
    else if(variable_state == SLEQP_ACTIVE_BOTH)
    {
      assert(sleqp_is_eq(ldval, udval, eps));

      SLEQP_CALL(sleqp_sparse_vector_push(rhs, variable_index, udval));
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
  const double zero_eps = sleqp_params_get_zero_eps(soc_data->params);

  SleqpSparseVec* lower_diff = soc_data->lower_diff;
  SleqpSparseVec* upper_diff = soc_data->upper_diff;

  SleqpSparseVec* rhs = soc_data->rhs;
  SleqpWorkingSet* working_set = sleqp_iterate_get_working_set(iterate);

  SLEQP_CALL(sleqp_sparse_vector_clear(soc_data->lower_diff));
  SLEQP_CALL(sleqp_sparse_vector_resize(soc_data->lower_diff, problem->num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_clear(soc_data->upper_diff));
  SLEQP_CALL(sleqp_sparse_vector_resize(soc_data->upper_diff, problem->num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(sleqp_iterate_get_cons_val(trial_iterate),
                                            problem->cons_lb,
                                            -1.,
                                            1.,
                                            zero_eps,
                                            lower_diff));

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(sleqp_iterate_get_cons_val(trial_iterate),
                                            problem->cons_ub,
                                            -1.,
                                            1.,
                                            zero_eps,
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

    const SLEQP_ACTIVE_STATE constraint_state = sleqp_working_set_get_constraint_state(working_set, i);

    if(constraint_state == SLEQP_INACTIVE)
    {
      continue;
    }

    const int constraint_index = sleqp_working_set_get_constraint_index(working_set, i);

    if(constraint_state == SLEQP_ACTIVE_UPPER && !sleqp_is_zero(udval, eps))
    {
      SLEQP_CALL(sleqp_sparse_vector_push(rhs, constraint_index, udval));
    }
    else if(constraint_state == SLEQP_ACTIVE_LOWER && !sleqp_is_zero(ldval, eps))
    {
      SLEQP_CALL(sleqp_sparse_vector_push(rhs, constraint_index, ldval));
    }
    else if(constraint_state == SLEQP_ACTIVE_BOTH)
    {
      assert(sleqp_is_eq(ldval, udval, eps));

      SLEQP_CALL(sleqp_sparse_vector_push(rhs, constraint_index, udval));
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
  SleqpWorkingSet* working_set = sleqp_iterate_get_working_set(iterate);

  const int working_set_size = sleqp_working_set_size(working_set);

  SLEQP_CALL(sleqp_sparse_vector_clear(rhs));

  SLEQP_CALL(sleqp_sparse_vector_resize(rhs, working_set_size));

  SLEQP_CALL(sleqp_sparse_vector_reserve(rhs, working_set_size));

  SLEQP_CALL(add_variable_entries(soc_data, iterate, trial_iterate));
  SLEQP_CALL(add_constraint_entries(soc_data, iterate, trial_iterate));

  SLEQP_CALL(sleqp_aug_jacobian_min_norm_solution(aug_jacobian, rhs, soc_direction));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_soc_data_free(SleqpSOCData** star)
{
  SleqpSOCData* soc_data = *star;

  if(!soc_data)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_sparse_vector_free(&soc_data->lower_diff));
  SLEQP_CALL(sleqp_sparse_vector_free(&soc_data->upper_diff));

  SLEQP_CALL(sleqp_sparse_vector_free(&soc_data->rhs));

  sleqp_free(star);

  return SLEQP_OKAY;
}
