#ifndef SLEQP_STEP_RULE_DIRECT_H
#define SLEQP_STEP_RULE_DIRECT_H

#include "step_rule.h"

/**
 * The default step rule, accepting iterates based on the reduction
 * ratio of the exact and model reduction.
 **/
SLEQP_RETCODE
sleqp_step_rule_direct_create(SleqpStepRule** star,
                              SleqpProblem* problem,
                              SleqpSettings* settings);

#endif /* SLEQP_STEP_RULE_DIRECT_H */
