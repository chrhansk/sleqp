#ifndef SLEQP_STEP_RULE_H
#define SLEQP_STEP_RULE_H

#include "iterate.h"
#include "options.h"
#include "types.h"

typedef struct SleqpStepRule SleqpStepRule;

typedef SLEQP_RETCODE (*SLEQP_STEP_RULE_APPLY)(double iterate_merit,
                                               double trial_exact_merit,
                                               double trial_model_merit,
                                               bool* accept_step,
                                               double* redution_ratio,
                                               void* step_data);

typedef SLEQP_RETCODE (*SLEQP_STEP_RULE_RESET)(void* step_data);

typedef SLEQP_RETCODE (*SLEQP_STEP_RULE_FREE)(void* step_data);

typedef struct
{
  SLEQP_STEP_RULE_APPLY rule_apply;
  SLEQP_STEP_RULE_RESET rule_reset;
  SLEQP_STEP_RULE_FREE rule_free;

} SleqpStepRuleCallbacks;

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_step_rule_create(SleqpStepRule** star,
                       SleqpProblem* problem,
                       SleqpParams* params,
                       SleqpStepRuleCallbacks* callbacks,
                       void* step_data);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_step_rule_apply(SleqpStepRule* rule,
                      double iterate_merit,
                      double trial_exact_merit,
                      double trial_model_merit,
                      bool* accept_step,
                      double* reduction_ratio);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_step_rule_reset(SleqpStepRule* rule);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_step_rule_capture(SleqpStepRule* rule);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_step_rule_release(SleqpStepRule** star);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_step_rule_create_default(SleqpStepRule** star,
                               SleqpProblem* problem,
                               SleqpParams* params,
                               const SleqpOptions* options);

#endif /* SLEQP_STEP_RULE_H */
