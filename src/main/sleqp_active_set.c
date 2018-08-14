#include "sleqp_active_set.h"

#include "sleqp_mem.h"

struct SleqpActiveSet
{
  SLEQP_ACTIVE_STATE* var_states;
  SLEQP_ACTIVE_STATE* cons_states;
};

SLEQP_RETCODE sleqp_active_set_create(SleqpActiveSet** star,
                                      SleqpProblem* problem)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpActiveSet* active_set = *star;

  SLEQP_CALL(sleqp_calloc(&active_set->var_states, problem->num_variables));
  SLEQP_CALL(sleqp_calloc(&active_set->cons_states, problem->num_constraints));

  return SLEQP_OKAY;
}

SLEQP_ACTIVE_STATE* sleqp_active_set_var_states(SleqpActiveSet* active_set)
{
  return active_set->var_states;
}

SLEQP_ACTIVE_STATE* sleqp_active_set_cons_states(SleqpActiveSet* active_set)
{
  return active_set->cons_states;
}

SLEQP_RETCODE sleqp_active_set_free(SleqpActiveSet** star)
{
  SleqpActiveSet* active_set = *star;

  sleqp_free(&active_set->cons_states);
  sleqp_free(&active_set->var_states);

  sleqp_free(star);

  return SLEQP_OKAY;
}
