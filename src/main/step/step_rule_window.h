#ifndef SLEQP_STEP_RULE_WINDOW_H
#define SLEQP_STEP_RULE_WINDOW_H

#include "step_rule.h"

/**
 * A non-monotone step rule, based on the maximum of the current
 * and a historic reduction ratio computed over a sliding window
 * with respect to the given window size.
 *
 * See "Trust-region methods", pp. 355
 **/
SLEQP_RETCODE
sleqp_step_rule_window_create(SleqpStepRule** star,
                              SleqpProblem* problem,
                              SleqpSettings* settings,
                              int window_size);

#endif /* SLEQP_STEP_RULE_WINDOW_H */
