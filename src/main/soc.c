#include "soc.h"

#include "cmp.h"
#include "fail.h"
#include "mem.h"
#include "settings.h"
#include "working_set.h"

#include "sparse/mat.h"

struct SleqpSOC
{
  int refcount;

  SleqpProblem* problem;
  SleqpSettings* settings;

  SleqpVec* upper_diff;
  SleqpVec* lower_diff;

  SleqpVec* soc_direction;
  SleqpVec* soc_corrected_direction;

  SleqpVec* rhs;
};

SLEQP_RETCODE
sleqp_soc_data_create(SleqpSOC** star,
                      SleqpProblem* problem,
                      SleqpSettings* settings)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpSOC* soc_data = *star;

  *soc_data          = (SleqpSOC){0};
  soc_data->refcount = 1;

  soc_data->problem = problem;
  SLEQP_CALL(sleqp_problem_capture(soc_data->problem));

  SLEQP_CALL(sleqp_settings_capture(settings));
  soc_data->settings = settings;

  const int num_variables   = sleqp_problem_num_vars(problem);
  const int num_constraints = sleqp_problem_num_cons(problem);

  SLEQP_CALL(sleqp_vec_create_empty(&soc_data->upper_diff, num_constraints));

  SLEQP_CALL(sleqp_vec_create_empty(&soc_data->lower_diff, num_constraints));

  SLEQP_CALL(sleqp_vec_create_empty(&soc_data->soc_direction, num_variables));

  SLEQP_CALL(
    sleqp_vec_create_empty(&soc_data->soc_corrected_direction, num_variables));

  SLEQP_CALL(sleqp_vec_create_empty(&soc_data->rhs, 0));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
add_variable_entries(SleqpSOC* soc_data,
                     const SleqpIterate* iterate,
                     const SleqpIterate* trial_iterate)
{
  SleqpProblem* problem = soc_data->problem;

  const int num_variables = sleqp_problem_num_vars(problem);

  const double eps = sleqp_settings_real_value(soc_data->settings, SLEQP_SETTINGS_REAL_EPS);

  const double zero_eps
    = sleqp_settings_real_value(soc_data->settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  SleqpVec* lower_diff = soc_data->lower_diff;
  SleqpVec* upper_diff = soc_data->upper_diff;

  SleqpVec* rhs                = soc_data->rhs;
  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);

  SLEQP_CALL(sleqp_vec_clear(soc_data->lower_diff));
  SLEQP_CALL(sleqp_vec_resize(soc_data->lower_diff, num_variables));

  SLEQP_CALL(sleqp_vec_clear(soc_data->upper_diff));
  SLEQP_CALL(sleqp_vec_resize(soc_data->upper_diff, num_variables));

  SLEQP_CALL(sleqp_vec_add_scaled(sleqp_iterate_primal(trial_iterate),
                                  sleqp_problem_vars_lb(problem),
                                  -1.,
                                  1.,
                                  zero_eps,
                                  lower_diff));

  SLEQP_CALL(sleqp_vec_add_scaled(sleqp_iterate_primal(trial_iterate),
                                  sleqp_problem_vars_ub(problem),
                                  -1.,
                                  1.,
                                  zero_eps,
                                  upper_diff));

  int k_ld = 0., k_ud = 0.;

  for (int i = 0; i < num_variables; ++i)
  {
    while (k_ld < lower_diff->nnz && lower_diff->indices[k_ld] < i)
    {
      ++k_ld;
    }

    while (k_ud < upper_diff->nnz && upper_diff->indices[k_ud] < i)
    {
      ++k_ud;
    }

    bool valid_ud = (k_ud < upper_diff->nnz && upper_diff->indices[k_ud] == i);
    bool valid_ld = (k_ld < lower_diff->nnz && lower_diff->indices[k_ld] == i);

    const double udval = valid_ud ? upper_diff->data[k_ud] : 0.;
    const double ldval = valid_ld ? lower_diff->data[k_ld] : 0.;

    const SLEQP_ACTIVE_STATE variable_state
      = sleqp_working_set_var_state(working_set, i);

    if (variable_state == SLEQP_INACTIVE)
    {
      continue;
    }

    const int variable_index = sleqp_working_set_var_index(working_set, i);

    if (variable_state == SLEQP_ACTIVE_UPPER && !sleqp_is_zero(udval, eps))
    {
      SLEQP_CALL(sleqp_vec_push(rhs, variable_index, udval));
    }
    else if (variable_state == SLEQP_ACTIVE_LOWER && !sleqp_is_zero(ldval, eps))
    {
      SLEQP_CALL(sleqp_vec_push(rhs, variable_index, ldval));
    }
    else if (variable_state == SLEQP_ACTIVE_BOTH)
    {
      sleqp_assert_is_eq(ldval, udval, eps);

      SLEQP_CALL(sleqp_vec_push(rhs, variable_index, udval));
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
add_constraint_entries(SleqpSOC* soc_data,
                       const SleqpIterate* iterate,
                       const SleqpIterate* trial_iterate)
{
  SleqpProblem* problem = soc_data->problem;

  const int num_constraints = sleqp_problem_num_cons(problem);

  const double eps = sleqp_settings_real_value(soc_data->settings, SLEQP_SETTINGS_REAL_EPS);

  const double zero_eps
    = sleqp_settings_real_value(soc_data->settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  SleqpVec* lower_diff = soc_data->lower_diff;
  SleqpVec* upper_diff = soc_data->upper_diff;

  SleqpVec* rhs                = soc_data->rhs;
  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);

  SLEQP_CALL(sleqp_vec_clear(soc_data->lower_diff));
  SLEQP_CALL(sleqp_vec_resize(soc_data->lower_diff, num_constraints));

  SLEQP_CALL(sleqp_vec_clear(soc_data->upper_diff));
  SLEQP_CALL(sleqp_vec_resize(soc_data->upper_diff, num_constraints));

  SLEQP_CALL(sleqp_vec_add_scaled(sleqp_iterate_cons_val(trial_iterate),
                                  sleqp_problem_cons_lb(problem),
                                  -1.,
                                  1.,
                                  zero_eps,
                                  lower_diff));

  SLEQP_CALL(sleqp_vec_add_scaled(sleqp_iterate_cons_val(trial_iterate),
                                  sleqp_problem_cons_ub(problem),
                                  -1.,
                                  1.,
                                  zero_eps,
                                  upper_diff));

  int k_ld = 0., k_ud = 0.;

  for (int i = 0; i < num_constraints; ++i)
  {
    while (k_ld < lower_diff->nnz && lower_diff->indices[k_ld] < i)
    {
      ++k_ld;
    }

    while (k_ud < upper_diff->nnz && upper_diff->indices[k_ud] < i)
    {
      ++k_ud;
    }

    bool valid_ud = (k_ud < upper_diff->nnz && upper_diff->indices[k_ud] == i);
    bool valid_ld = (k_ld < lower_diff->nnz && lower_diff->indices[k_ld] == i);

    const double udval = valid_ud ? upper_diff->data[k_ud] : 0.;
    const double ldval = valid_ld ? lower_diff->data[k_ld] : 0.;

    const SLEQP_ACTIVE_STATE constraint_state
      = sleqp_working_set_cons_state(working_set, i);

    if (constraint_state == SLEQP_INACTIVE)
    {
      continue;
    }

    const int constraint_index = sleqp_working_set_cons_index(working_set, i);

    if (constraint_state == SLEQP_ACTIVE_UPPER && !sleqp_is_zero(udval, eps))
    {
      SLEQP_CALL(sleqp_vec_push(rhs, constraint_index, udval));
    }
    else if (constraint_state == SLEQP_ACTIVE_LOWER
             && !sleqp_is_zero(ldval, eps))
    {
      SLEQP_CALL(sleqp_vec_push(rhs, constraint_index, ldval));
    }
    else if (constraint_state == SLEQP_ACTIVE_BOTH)
    {
      sleqp_assert_is_eq(ldval, udval, eps);

      SLEQP_CALL(sleqp_vec_push(rhs, constraint_index, udval));
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_soc_compute_correction(SleqpSOC* soc_data,
                             SleqpAugJac* aug_jac,
                             const SleqpIterate* iterate,
                             const SleqpIterate* trial_iterate,
                             SleqpVec* soc_direction)
{
  SleqpVec* rhs                = soc_data->rhs;
  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);

  const int working_set_size = sleqp_working_set_size(working_set);

  SLEQP_CALL(sleqp_vec_clear(rhs));

  SLEQP_CALL(sleqp_vec_resize(rhs, working_set_size));

  SLEQP_CALL(sleqp_vec_reserve(rhs, working_set_size));

  SLEQP_CALL(add_variable_entries(soc_data, iterate, trial_iterate));
  SLEQP_CALL(add_constraint_entries(soc_data, iterate, trial_iterate));

  SLEQP_CALL(sleqp_aug_jac_solve_min_norm(aug_jac, rhs, soc_direction));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_soc_compute_step(SleqpSOC* soc_data,
                       SleqpAugJac* aug_jac,
                       const SleqpIterate* iterate,
                       const SleqpVec* trial_step,
                       const SleqpIterate* trial_iterate,
                       SleqpVec* soc_step)
{
  assert(trial_step != soc_step);

  SleqpProblem* problem = soc_data->problem;

  SleqpVec* trial_point = sleqp_iterate_primal(trial_iterate);

  const double zero_eps
    = sleqp_settings_real_value(soc_data->settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  double max_step_length = 1.;

  SLEQP_CALL(sleqp_soc_compute_correction(soc_data,
                                          aug_jac,
                                          iterate,
                                          trial_iterate,
                                          soc_data->soc_direction));

  SLEQP_CALL(sleqp_max_step_length(trial_point,
                                   soc_data->soc_direction,
                                   sleqp_problem_vars_lb(problem),
                                   sleqp_problem_vars_ub(problem),
                                   &max_step_length));

  SLEQP_CALL(sleqp_vec_add_scaled(trial_step,
                                  soc_data->soc_direction,
                                  1.,
                                  max_step_length,
                                  zero_eps,
                                  soc_step));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_soc_compute_trial_point(SleqpSOC* soc_data,
                              SleqpAugJac* aug_jac,
                              const SleqpIterate* iterate,
                              const SleqpVec* trial_step,
                              const SleqpIterate* trial_iterate,
                              SleqpVec* soc_trial_point,
                              double* soc_step_norm)
{
  SleqpProblem* problem = soc_data->problem;

  SleqpVec* current_point = sleqp_iterate_primal(iterate);
  SleqpVec* trial_point   = sleqp_iterate_primal(trial_iterate);

  const double zero_eps
    = sleqp_settings_real_value(soc_data->settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  SLEQP_CALL(sleqp_soc_compute_correction(soc_data,
                                          aug_jac,
                                          iterate,
                                          trial_iterate,
                                          soc_data->soc_direction));

  double max_step_length = 1.;

  SLEQP_CALL(sleqp_max_step_length(trial_point,
                                   soc_data->soc_direction,
                                   sleqp_problem_vars_lb(problem),
                                   sleqp_problem_vars_ub(problem),
                                   &max_step_length));

  SLEQP_CALL(sleqp_vec_add_scaled(trial_step,
                                  soc_data->soc_direction,
                                  1.,
                                  max_step_length,
                                  zero_eps,
                                  soc_data->soc_corrected_direction));

  (*soc_step_norm) = sleqp_vec_norm(soc_data->soc_corrected_direction);

  SLEQP_CALL(sleqp_vec_add_scaled(current_point,
                                  soc_data->soc_corrected_direction,
                                  1.,
                                  1.,
                                  zero_eps,
                                  soc_trial_point));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
soc_data_free(SleqpSOC** star)
{
  SleqpSOC* soc_data = *star;

  if (!soc_data)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_vec_free(&soc_data->rhs));

  SLEQP_CALL(sleqp_vec_free(&soc_data->soc_corrected_direction));
  SLEQP_CALL(sleqp_vec_free(&soc_data->soc_direction));

  SLEQP_CALL(sleqp_vec_free(&soc_data->lower_diff));
  SLEQP_CALL(sleqp_vec_free(&soc_data->upper_diff));

  SLEQP_CALL(sleqp_settings_release(&soc_data->settings));

  SLEQP_CALL(sleqp_problem_release(&soc_data->problem));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_soc_data_capture(SleqpSOC* soc_data)
{
  ++soc_data->refcount;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_soc_data_release(SleqpSOC** star)
{
  SleqpSOC* soc_data = *star;

  if (!soc_data)
  {
    return SLEQP_OKAY;
  }

  if (--soc_data->refcount == 0)
  {
    SLEQP_CALL(soc_data_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
