#include "step_rule_direct.h"

#include "mem.h"
#include "util.h"

typedef struct
{
  SleqpSettings* settings;
} StepRule;

static SLEQP_RETCODE
step_rule_direct_apply(double iterate_merit,
                       double trial_exact_merit,
                       double trial_model_merit,
                       bool* accept_step,
                       double* reduction_ratio,
                       void* step_data)
{
  StepRule* step_rule = (StepRule*)step_data;

  const double exact_reduction = iterate_merit - trial_exact_merit;
  const double model_reduction = iterate_merit - trial_model_merit;

  *reduction_ratio = sleqp_reduction_ratio(exact_reduction, model_reduction);

  const double accepted_reduction
    = sleqp_settings_real_value(step_rule->settings, SLEQP_SETTINGS_REAL_ACCEPTED_REDUCTION);

  *accept_step = (*reduction_ratio >= accepted_reduction);

  sleqp_log_debug("Step with reduction ratio %e, accepted: %s",
                  *reduction_ratio,
                  sleqp_bool_string(*accept_step));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
step_rule_direct_free(void* step_data)
{
  StepRule* step_rule = (StepRule*)step_data;

  SLEQP_CALL(sleqp_settings_release(&step_rule->settings));

  sleqp_free(&step_data);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_step_rule_direct_create(SleqpStepRule** star,
                              SleqpProblem* problem,
                              SleqpSettings* settings)
{
  SleqpStepRuleCallbacks callbacks = {.rule_apply = step_rule_direct_apply,
                                      .rule_reset = NULL,
                                      .rule_free  = step_rule_direct_free};

  StepRule* step_rule = NULL;

  SLEQP_CALL(sleqp_malloc(&step_rule));

  SLEQP_CALL(sleqp_settings_capture(settings));
  step_rule->settings = settings;

  SLEQP_CALL(sleqp_step_rule_create(star,
                                    problem,
                                    settings,
                                    &callbacks,
                                    (void*)step_rule));

  return SLEQP_OKAY;
}
