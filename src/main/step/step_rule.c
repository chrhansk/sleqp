#include "step_rule.h"

#include <assert.h>
#include <float.h>

#include "step/step_rule_direct.h"
#include "step/step_rule_minstep.h"
#include "step/step_rule_window.h"

#include "cmp.h"
#include "mem.h"

static const int window_size = 25;
static const int step_count = 2;

static const double eps_factor = 10.;

struct SleqpStepRule
{
  int refcount;

  SleqpStepRuleCallbacks callbacks;
  SleqpProblem* problem;

  void* step_data;
};

SLEQP_RETCODE sleqp_step_rule_create(SleqpStepRule** star,
                                     SleqpProblem* problem,
                                     SleqpParams* params,
                                     SleqpStepRuleCallbacks* callbacks,
                                     void* step_data)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpStepRule* rule = *star;

  *rule = (SleqpStepRule) {0};

  rule->refcount = 1;

  rule->callbacks = *callbacks;

  SLEQP_CALL(sleqp_problem_capture(problem));

  rule->problem = problem;
  rule->step_data = step_data;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_step_rule_apply(SleqpStepRule* rule,
                                    double iterate_merit,
                                    double trial_exact_merit,
                                    double trial_model_merit,
                                    bool* accept_step,
                                    double* redution_ratio)
{
  SLEQP_CALL(rule->callbacks.rule_apply(iterate_merit,
                                        trial_exact_merit,
                                        trial_model_merit,
                                        accept_step,
                                        redution_ratio,
                                        rule->step_data));

  return SLEQP_OKAY;
}

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_step_rule_reset(SleqpStepRule* rule)
{
  if(rule->callbacks.rule_reset)
  {
    SLEQP_CALL(rule->callbacks.rule_reset(rule->step_data));
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_step_rule_capture(SleqpStepRule* rule)
{
  ++rule->refcount;

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE step_rule_free(SleqpStepRule** star)
{
  SleqpStepRule* rule = *star;

  if(!rule)
  {
    return SLEQP_OKAY;
  }

  if(rule->callbacks.rule_free)
  {
    SLEQP_CALL(rule->callbacks.rule_free(rule->step_data));
  }

  rule->step_data = NULL;

  SLEQP_CALL(sleqp_problem_release(&rule->problem));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_step_rule_release(SleqpStepRule** star)
{
  SleqpStepRule* rule = *star;

  if(!rule)
  {
    return SLEQP_OKAY;
  }

  if(--rule->refcount == 0)
  {
    SLEQP_CALL(step_rule_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_step_rule_create_default(SleqpStepRule** star,
                                             SleqpProblem* problem,
                                             SleqpParams* params,
                                             const SleqpOptions* options)
{
  SLEQP_STEP_RULE step_rule = sleqp_options_get_int(options,
                                                    SLEQP_OPTION_INT_STEP_RULE);

  if(step_rule == SLEQP_STEP_RULE_DIRECT)
  {
    SLEQP_CALL(sleqp_step_rule_direct_create(star,
                                             problem,
                                             params));
  }
  else if(step_rule == SLEQP_STEP_RULE_WINDOW)
  {
    SLEQP_CALL(sleqp_step_rule_window_create(star,
                                             problem,
                                             params,
                                             window_size));
  }
  else
  {
    assert(step_rule == SLEQP_STEP_RULE_MINSTEP);

    SLEQP_CALL(sleqp_step_rule_minstep_create(star,
                                              problem,
                                              params,
                                              step_count));
  }

  return SLEQP_OKAY;
}

double sleqp_step_rule_reduction_ratio(const double exact_reduction,
                                       const double model_reduction)
{
  const double eps = eps_factor * DBL_EPSILON;

  const double corr_model_reduction = model_reduction - eps;
  const double corr_exact_reduction = exact_reduction - eps;

  // Safeguard against roundoff errors
  if(SLEQP_ABS(corr_model_reduction) <= eps &&
     SLEQP_ABS(corr_exact_reduction) <= eps)
  {
    return 1.;
  }

  return corr_exact_reduction / corr_model_reduction;
}
