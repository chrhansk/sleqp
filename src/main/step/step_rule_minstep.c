#include "step_rule_minstep.h"

#include <assert.h>

#include "cmp.h"
#include "log.h"
#include "mem.h"

typedef struct
{
  SleqpParams* params;

  bool init;

  // min merit ever seen
  double min_exact_merit;

  // merit at reference iteration
  double ref_exact_merit;

  // max merit seen since reaching
  // current reference iteration
  double max_exact_merit;

  // sum of model decreases of accepted steps
  // since the model iterate
  double model_decrease_sum_ref;

  // sum of model decreases of accepted steps
  // since the max iterate
  double model_decrease_sum_max;

  double last_exact_merit;

  int step_count;
  int max_step_count;

} StepRule;

static SLEQP_RETCODE
step_rule_init(StepRule* step_rule, double current_merit)
{
  assert(!step_rule->init);

  step_rule->step_count      = 0;
  step_rule->min_exact_merit = current_merit;
  step_rule->ref_exact_merit = current_merit;
  step_rule->max_exact_merit = current_merit;

  step_rule->model_decrease_sum_ref = 0.;
  step_rule->model_decrease_sum_max = 0.;

  step_rule->init = true;

  return SLEQP_OKAY;
}

static double
compute_historic_reduction_ratio(const StepRule* step_rule,
                                 const double trial_exact_merit,
                                 const double model_reduction)
{
  const double total_model_reduction
    = step_rule->model_decrease_sum_ref + model_reduction;

  return (step_rule->ref_exact_merit - trial_exact_merit)
         / total_model_reduction;
}

static SLEQP_RETCODE
step_rule_minstep_apply(double iterate_merit,
                        double trial_exact_merit,
                        double trial_model_merit,
                        bool* accept_step,
                        double* reduction_ratio,
                        void* step_data)
{
  StepRule* step_rule = (StepRule*)step_data;

  const double exact_reduction = iterate_merit - trial_exact_merit;
  const double model_reduction = iterate_merit - trial_model_merit;

  step_rule->last_exact_merit = trial_exact_merit;

  if (!step_rule->init)
  {
    SLEQP_CALL(step_rule_init(step_rule, iterate_merit));
  }
  else
  {
    assert(step_rule->max_exact_merit >= step_rule->min_exact_merit);
  }

  const double historic_reduction_ratio
    = compute_historic_reduction_ratio(step_rule,
                                       trial_exact_merit,
                                       model_reduction);

  double current_reduction_ratio = 1.;

  if (exact_reduction != model_reduction)
  {
    current_reduction_ratio = exact_reduction / model_reduction;
  }

  *reduction_ratio
    = SLEQP_MAX(current_reduction_ratio, historic_reduction_ratio);

  const double accepted_reduction
    = sleqp_params_get(step_rule->params, SLEQP_PARAM_ACCEPTED_REDUCTION);

  *accept_step = (*reduction_ratio >= accepted_reduction);

  if (!(*accept_step))
  {
    return SLEQP_OKAY;
  }

  if (current_reduction_ratio < accepted_reduction)
  {
    sleqp_log_debug("Accepting due to historic reduction ratio");
  }

  if (step_rule->last_exact_merit < trial_exact_merit)
  {
    sleqp_log_debug("Accepting an iterate with increasing merit");
  }

  step_rule->model_decrease_sum_ref += model_reduction;
  step_rule->model_decrease_sum_max += model_reduction;

  // found a new minimum
  if (iterate_merit < step_rule->min_exact_merit)
  {
    sleqp_log_debug("Found new minimum with merit %e", iterate_merit);

    step_rule->min_exact_merit = iterate_merit;
    step_rule->max_exact_merit = iterate_merit;

    step_rule->model_decrease_sum_ref = 0.;
    step_rule->model_decrease_sum_max = 0.;
    step_rule->step_count             = 0;
  }
  else
  {
    ++(step_rule->step_count);
  }

  if (iterate_merit > step_rule->max_exact_merit)
  {
    sleqp_log_debug("Found new local maximum with merit %e", iterate_merit);

    step_rule->max_exact_merit        = iterate_merit;
    step_rule->model_decrease_sum_max = 0.;
  }

  if (step_rule->step_count == step_rule->max_step_count)
  {
    sleqp_log_debug("Reached step count limit, resetting reference merit to %e",
                    step_rule->max_exact_merit);

    step_rule->ref_exact_merit        = step_rule->max_exact_merit;
    step_rule->model_decrease_sum_ref = step_rule->model_decrease_sum_max;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
step_rule_minstep_reset(void* step_data)
{
  StepRule* step_rule = (StepRule*)step_data;

  step_rule->init = false;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
step_rule_minstep_free(void* step_data)
{
  StepRule* step_rule = (StepRule*)step_data;

  SLEQP_CALL(sleqp_params_release(&step_rule->params));

  sleqp_free(&step_rule);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_step_rule_minstep_create(SleqpStepRule** star,
                               SleqpProblem* problem,
                               SleqpParams* params,
                               int step_count)
{
  assert(step_count > 0);

  SleqpStepRuleCallbacks callbacks = {.rule_apply = step_rule_minstep_apply,
                                      .rule_reset = step_rule_minstep_reset,
                                      .rule_free  = step_rule_minstep_free};

  StepRule* step_rule = NULL;

  SLEQP_CALL(sleqp_malloc(&step_rule));

  SLEQP_CALL(sleqp_params_capture(params));
  step_rule->params = params;

  step_rule->init           = false;
  step_rule->max_step_count = step_count;

  SLEQP_CALL(sleqp_step_rule_create(star,
                                    problem,
                                    params,
                                    &callbacks,
                                    (void*)step_rule));

  return SLEQP_OKAY;
}
