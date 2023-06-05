#include "polish.h"

#include "cmp.h"
#include "fail.h"
#include "log.h"
#include "mem.h"
#include "settings.h"
#include "working_set.h"

struct SleqpPolishing
{
  int refcount;

  SleqpProblem* problem;
  SleqpSettings* settings;
  SleqpWorkingSet* working_set;
};

SLEQP_RETCODE
sleqp_polishing_create(SleqpPolishing** star,
                       SleqpProblem* problem,
                       SleqpSettings* settings)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpPolishing* polishing = *star;

  *polishing = (SleqpPolishing){0};

  polishing->refcount = 1;

  SLEQP_CALL(sleqp_problem_capture(problem));
  polishing->problem = problem;

  SLEQP_CALL(sleqp_settings_capture(settings));
  polishing->settings = settings;

  SLEQP_CALL(sleqp_working_set_create(&polishing->working_set, problem));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
polish_inactive_range(SleqpPolishing* polishing,
                      const SleqpWorkingSet* source_set,
                      const SleqpVec* lb,
                      const SleqpVec* primal,
                      const SleqpVec* ub,
                      const SleqpVec* dual,
                      bool constraints,
                      int* num_removed)
{
  SleqpSettings* settings = polishing->settings;

  SleqpWorkingSet* target_set = polishing->working_set;

  (*num_removed) = 0;

  const double feas_eps = sleqp_settings_real_value(settings, SLEQP_SETTINGS_REAL_FEAS_TOL);

  const int size = dual->dim;
  int k_d        = 0;
  int k_lb = 0, k_ub = 0, k_p = 0;

  for (int j = 0; j < size; ++j)
  {
    SLEQP_ACTIVE_STATE state
      = sleqp_working_set_state(source_set, constraints, j);

    if (state == SLEQP_INACTIVE)
    {
      continue;
    }

    while (k_d < dual->nnz && dual->indices[k_d] < j)
    {
      ++k_d;
    }

    if (k_d < dual->nnz && dual->indices[k_d] == j)
    {
      SLEQP_CALL(sleqp_working_set_add(target_set, j, constraints, state));

      continue;
    }

    while (k_lb < lb->nnz && lb->indices[k_lb] < j)
    {
      ++k_lb;
    }

    while (k_ub < ub->nnz && ub->indices[k_ub] < j)
    {
      ++k_ub;
    }

    while (k_p < primal->nnz && primal->indices[k_p] < j)
    {
      ++k_p;
    }

    const bool valid_lb = (k_lb < lb->nnz && lb->indices[k_lb] == j);
    const double lb_val = valid_lb ? lb->data[k_lb] : 0.;

    const bool valid_ub = (k_ub < ub->nnz && ub->indices[k_ub] == j);
    const double ub_val = valid_ub ? ub->data[k_ub] : 0.;

    const bool valid_p = (k_p < primal->nnz && primal->indices[k_p] == j);
    const double p_val = valid_p ? primal->data[k_p] : 0.;

    const bool active_upper = sleqp_is_eq(p_val, ub_val, feas_eps);
    const bool active_lower = sleqp_is_eq(lb_val, p_val, feas_eps);

    bool inactive = !(active_lower || active_upper);

    if (inactive)
    {
      ++(*num_removed);
    }
    else
    {
      SLEQP_CALL(sleqp_working_set_add(target_set, j, constraints, state));
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
polish_zero_dual_range(SleqpPolishing* polishing,
                       const SleqpWorkingSet* source_set,
                       const SleqpVec* dual,
                       bool constraints,
                       int* num_removed)
{
  SleqpWorkingSet* target_set = polishing->working_set;

  (*num_removed) = 0;

  const int size = dual->dim;

  int k = 0;

  for (int i = 0; i < size; ++i)
  {
    SLEQP_ACTIVE_STATE state
      = sleqp_working_set_state(source_set, constraints, i);

    if (state == SLEQP_INACTIVE)
    {
      continue;
    }

    while (k < dual->nnz && dual->indices[k] < i)
    {
      ++k;
    }

    if (k < dual->nnz && dual->indices[k] == i)
    {
      SLEQP_CALL(sleqp_working_set_add(target_set, i, constraints, state));
    }
    else
    {
      ++(*num_removed);
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
polish_inactive(SleqpPolishing* polishing, SleqpIterate* iterate)
{
  SleqpProblem* problem = polishing->problem;

  int num_removed_vars = 0, num_removed_cons = 0;

  SleqpWorkingSet* source_set = sleqp_iterate_working_set(iterate);
  SleqpWorkingSet* target_set = polishing->working_set;

  SLEQP_CALL(sleqp_working_set_reset(target_set));

  SLEQP_CALL(polish_inactive_range(polishing,
                                   source_set,
                                   sleqp_problem_vars_lb(problem),
                                   sleqp_iterate_primal(iterate),
                                   sleqp_problem_vars_ub(problem),
                                   sleqp_iterate_vars_dual(iterate),
                                   false,
                                   &num_removed_vars));

  SLEQP_CALL(polish_inactive_range(polishing,
                                   source_set,
                                   sleqp_problem_cons_lb(problem),
                                   sleqp_iterate_cons_val(iterate),
                                   sleqp_problem_cons_ub(problem),
                                   sleqp_iterate_cons_dual(iterate),
                                   true,
                                   &num_removed_cons));

  sleqp_log_debug("Polishing removed %d variables and %d constraints",
                  num_removed_vars,
                  num_removed_cons);

  SLEQP_CALL(sleqp_working_set_copy(target_set, source_set));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
polish_zero_dual(SleqpPolishing* polishing, SleqpIterate* iterate)
{
  int num_removed_vars = 0, num_removed_cons = 0;

  SleqpWorkingSet* source_set = sleqp_iterate_working_set(iterate);
  SleqpWorkingSet* target_set = polishing->working_set;

  SLEQP_CALL(sleqp_working_set_reset(target_set));

  SLEQP_CALL(polish_zero_dual_range(polishing,
                                    source_set,
                                    sleqp_iterate_vars_dual(iterate),
                                    false,
                                    &num_removed_vars));

  SLEQP_CALL(polish_zero_dual_range(polishing,
                                    source_set,
                                    sleqp_iterate_cons_dual(iterate),
                                    true,
                                    &num_removed_cons));

  sleqp_log_debug("Polishing removed %d variables and %d constraints",
                  num_removed_vars,
                  num_removed_cons);

  SLEQP_CALL(sleqp_working_set_copy(target_set, source_set));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_polishing_polish(SleqpPolishing* polishing,
                       SleqpIterate* iterate,
                       SLEQP_POLISHING_TYPE polishing_type)
{
  if (polishing_type == SLEQP_POLISHING_NONE)
  {
    return SLEQP_OKAY;
  }
  else if (polishing_type == SLEQP_POLISHING_ZERO_DUAL)
  {
    SLEQP_CALL(polish_zero_dual(polishing, iterate));
  }
  else
  {
    assert(polishing_type == SLEQP_POLISHING_INACTIVE);

    SLEQP_CALL(polish_inactive(polishing, iterate));
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_polishing_capture(SleqpPolishing* polishing)
{
  ++polishing->refcount;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
polishing_free(SleqpPolishing** star)
{
  SleqpPolishing* polishing = *star;

  SLEQP_CALL(sleqp_working_set_release(&polishing->working_set));

  SLEQP_CALL(sleqp_settings_release(&polishing->settings));

  SLEQP_CALL(sleqp_problem_release(&polishing->problem));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_polishing_release(SleqpPolishing** star)
{
  SleqpPolishing* polishing = *star;

  if (!polishing)
  {
    return SLEQP_OKAY;
  }

  if (--polishing->refcount == 0)
  {
    SLEQP_CALL(polishing_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
