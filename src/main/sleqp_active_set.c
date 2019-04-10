#include "sleqp_active_set.h"

#include "sleqp_mem.h"

#define INACTIVE_INDEX -1

struct SleqpActiveSet
{
  SleqpProblem* problem;

  SLEQP_ACTIVE_STATE* variable_states;
  SLEQP_ACTIVE_STATE* constraint_states;

  int num_variables;
  int num_constraints;

  int max_set_size;

  int num_active_constraints;
  int num_active_variables;

  // a mapping of 0..num_variables - 1 -> pos in the active set or -1
  int* variable_indices;

  // a mapping of 0..num_constraints - 1 -> pos in the active set or -1
  int* constraint_indices;

  int* content_indices;

};

SLEQP_RETCODE sleqp_active_set_create(SleqpActiveSet** star,
                                      SleqpProblem* problem)
{
  SLEQP_CALL(sleqp_malloc(star));

  int num_variables = problem->num_variables;
  int num_constraints = problem->num_constraints;

  SleqpActiveSet* active_set = *star;

  active_set->problem = problem;

  SLEQP_CALL(sleqp_calloc(&active_set->variable_states, num_variables));
  SLEQP_CALL(sleqp_calloc(&active_set->constraint_states, num_constraints));

  active_set->num_variables = num_variables;
  active_set->num_constraints = num_constraints;

  active_set->max_set_size = num_variables + num_constraints;

  SLEQP_CALL(sleqp_calloc(&active_set->variable_indices, num_variables));
  SLEQP_CALL(sleqp_calloc(&active_set->constraint_indices, num_constraints));

  SLEQP_CALL(sleqp_calloc(&active_set->content_indices,
                          active_set->max_set_size));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_active_set_reset(SleqpActiveSet* active_set)
{
  for(int j = 0; j < active_set->num_variables; ++j)
  {
    active_set->variable_states[j] = SLEQP_INACTIVE;
    active_set->variable_indices[j] = -1;
  }

  for(int i = 0; i < active_set->num_constraints; ++i)
  {
    active_set->constraint_states[i] = SLEQP_INACTIVE;
    active_set->constraint_indices[i] = -1;
  }

  active_set->num_active_constraints = 0;
  active_set->num_active_variables = 0;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_active_set_add_variable(SleqpActiveSet* active_set,
                                            int index,
                                            SLEQP_ACTIVE_STATE state)
{
  assert(active_set->num_active_constraints == 0);

  assert(sleqp_active_set_get_variable_state(active_set, index) == SLEQP_INACTIVE);
  assert(state != SLEQP_INACTIVE);

  const int size = sleqp_active_set_size(active_set);

  active_set->variable_indices[index] = (active_set->num_active_variables)++;

  active_set->content_indices[size] = index;

  active_set->variable_states[index] = state;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_active_set_add_constraint(SleqpActiveSet* active_set,
                                              int index,
                                              SLEQP_ACTIVE_STATE state)
{
  assert(sleqp_active_set_get_constraint_state(active_set, index) == SLEQP_INACTIVE);
  assert(state != SLEQP_INACTIVE);

  const int size = sleqp_active_set_size(active_set);

  active_set->constraint_indices[index] = active_set->num_active_variables
    + (active_set->num_active_constraints)++;

  active_set->content_indices[size] = active_set->num_variables + index;

  active_set->constraint_states[index] = state;

  return SLEQP_OKAY;
}

int sleqp_active_set_get_constraint_index(SleqpActiveSet* active_set,
                                          int index)
{
  assert(index < active_set->num_constraints);

  return active_set->constraint_indices[index];
}

int sleqp_active_set_get_variable_index(SleqpActiveSet* active_set,
                                        int index)
{
  assert(index < active_set->num_variables);

  return active_set->variable_indices[index];
}

int sleqp_active_set_get_content(SleqpActiveSet* active_set,
                                 int index)
{
  return active_set->content_indices[index];
}

SLEQP_ACTIVE_STATE* sleqp_active_set_variable_states(SleqpActiveSet* active_set)
{
  return active_set->variable_states;
}

SLEQP_ACTIVE_STATE* sleqp_active_set_constraint_states(SleqpActiveSet* active_set)
{
  return active_set->constraint_states;
}

SLEQP_ACTIVE_STATE sleqp_active_set_get_variable_state(SleqpActiveSet* active_set,
                                                       int j)
{
  assert(j < active_set->num_variables);
  return active_set->variable_states[j];
}

SLEQP_ACTIVE_STATE sleqp_active_set_get_constraint_state(SleqpActiveSet* active_set,
                                                         int i)
{
  assert(i < active_set->num_constraints);
  return active_set->constraint_states[i];
}

int sleqp_active_set_num_active_vars(SleqpActiveSet* active_set)
{
  return active_set->num_active_variables;
}

int sleqp_active_set_num_active_conss(SleqpActiveSet* active_set)
{
  return active_set->num_active_constraints;
}

int sleqp_active_set_size(SleqpActiveSet* active_set)
{
  return active_set->num_active_constraints + active_set->num_active_variables;
}

SLEQP_RETCODE sleqp_active_set_fprintf(SleqpActiveSet* active_set,
                                       FILE* output)
{
  SleqpProblem* problem = active_set->problem;

  int num_variables = problem->num_variables;
  int num_constraints = problem->num_constraints;


  fprintf(output,
          "Active set, variables: %d, constraints: %d\n",
          num_variables,
          num_constraints);

  const char* state_names[] = {[SLEQP_INACTIVE] = "inactive",
                               [SLEQP_ACTIVE_UPPER] = "upper",
                               [SLEQP_ACTIVE_LOWER] = "lower",
                               [SLEQP_ACTIVE_BOTH] = "active"};

  for(int j = 0; j < num_variables; ++j)
  {
    SLEQP_ACTIVE_STATE state = sleqp_active_set_get_variable_state(active_set, j);

    fprintf(output,
            "State of variable %d: %s\n",
            j,
            state_names[state]);
  }

  for(int i = 0; i < num_constraints; ++i)
  {
    SLEQP_ACTIVE_STATE state = sleqp_active_set_get_constraint_state(active_set, i);

    fprintf(output,
            "State of constraint %d: %s\n",
            i,
            state_names[state]);
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_active_set_copy(SleqpActiveSet* source,
                                    SleqpActiveSet* target)
{
  const int num_variables = source->problem->num_variables;
  const int num_constraints = source->problem->num_constraints;
  const int max_set_size = source->max_set_size;

  assert(num_variables == target->problem->num_variables);
  assert(num_constraints == target->problem->num_constraints);
  assert(max_set_size == target->max_set_size);

  for(int j = 0; j < num_variables; ++j)
  {
    target->variable_states[j] = source->variable_states[j];
  }

  for(int i = 0; i < num_constraints; ++i)
  {
    target->constraint_states[i] = source->constraint_states[i];
  }

  target->num_variables = source->num_variables;
  target->num_constraints = source->num_constraints;

  target->max_set_size = source->max_set_size;

  target->num_active_constraints = source->num_active_constraints;
  target->num_active_variables = source->num_active_variables;

  for(int j = 0; j < num_variables; ++j)
  {
    target->variable_indices[j] = source->variable_indices[j];
  }

  for(int i = 0; i < num_constraints; ++i)
  {
    target->constraint_indices[i] = source->constraint_indices[i];
  }

  for(int k = 0; k < max_set_size; ++k)
  {
    target->content_indices[k] = source->content_indices[k];
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_active_set_free(SleqpActiveSet** star)
{
  SleqpActiveSet* active_set = *star;

  if(!active_set)
  {
    return SLEQP_OKAY;
  }

  sleqp_free(&active_set->content_indices);

  sleqp_free(&active_set->constraint_indices);
  sleqp_free(&active_set->variable_indices);

  sleqp_free(&active_set->constraint_states);
  sleqp_free(&active_set->variable_states);

  sleqp_free(star);

  return SLEQP_OKAY;
}
