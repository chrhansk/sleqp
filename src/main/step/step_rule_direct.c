#include "step_rule_direct.h"

#include "mem.h"

typedef struct
{
  SleqpParams* params;
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

  *reduction_ratio
    = sleqp_step_rule_reduction_ratio(exact_reduction, model_reduction);

  const double accepted_reduction
    = sleqp_params_get(step_rule->params, SLEQP_PARAM_ACCEPTED_REDUCTION);

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

  SLEQP_CALL(sleqp_params_release(&step_rule->params));

  sleqp_free(&step_data);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_step_rule_direct_create(SleqpStepRule** star,
                              SleqpProblem* problem,
                              SleqpParams* params)
{
  SleqpStepRuleCallbacks callbacks = {.rule_apply = step_rule_direct_apply,
                                      .rule_reset = NULL,
                                      .rule_free  = step_rule_direct_free};

  StepRule* step_rule = NULL;

  SLEQP_CALL(sleqp_malloc(&step_rule));

  SLEQP_CALL(sleqp_params_capture(params));
  step_rule->params = params;

  SLEQP_CALL(sleqp_step_rule_create(star,
                                    problem,
                                    params,
                                    &callbacks,
                                    (void*)step_rule));

  return SLEQP_OKAY;
}
