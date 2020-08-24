#include "sleqp_working_set.h"

#include "sleqp_mem.h"

struct SleqpWorkingSet
{
  int refcount;

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

SLEQP_RETCODE sleqp_working_set_create(SleqpWorkingSet** star,
                                      SleqpProblem* problem)
{
  SLEQP_CALL(sleqp_malloc(star));

  int num_variables = problem->num_variables;
  int num_constraints = problem->num_constraints;

  SleqpWorkingSet* working_set = *star;

  working_set->refcount = 1;
  working_set->problem = problem;

  SLEQP_CALL(sleqp_calloc(&working_set->variable_states, num_variables));
  SLEQP_CALL(sleqp_calloc(&working_set->constraint_states, num_constraints));

  working_set->num_variables = num_variables;
  working_set->num_constraints = num_constraints;

  working_set->max_set_size = num_variables + num_constraints;

  SLEQP_CALL(sleqp_calloc(&working_set->variable_indices, num_variables));
  SLEQP_CALL(sleqp_calloc(&working_set->constraint_indices, num_constraints));

  SLEQP_CALL(sleqp_calloc(&working_set->content_indices,
                          working_set->max_set_size));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_working_set_reset(SleqpWorkingSet* working_set)
{
  for(int j = 0; j < working_set->num_variables; ++j)
  {
    working_set->variable_states[j] = SLEQP_INACTIVE;
    working_set->variable_indices[j] = -1;
  }

  for(int i = 0; i < working_set->num_constraints; ++i)
  {
    working_set->constraint_states[i] = SLEQP_INACTIVE;
    working_set->constraint_indices[i] = -1;
  }

  working_set->num_active_constraints = 0;
  working_set->num_active_variables = 0;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_working_set_add_variable(SleqpWorkingSet* working_set,
                                            int index,
                                            SLEQP_ACTIVE_STATE state)
{
  assert(working_set->num_active_constraints == 0);

  assert(sleqp_working_set_get_variable_state(working_set, index) == SLEQP_INACTIVE);
  assert(state != SLEQP_INACTIVE);

  const int size = sleqp_working_set_size(working_set);

  working_set->variable_indices[index] = (working_set->num_active_variables)++;

  working_set->content_indices[size] = index;

  working_set->variable_states[index] = state;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_working_set_add_constraint(SleqpWorkingSet* working_set,
                                              int index,
                                              SLEQP_ACTIVE_STATE state)
{
  assert(sleqp_working_set_get_constraint_state(working_set, index) == SLEQP_INACTIVE);
  assert(state != SLEQP_INACTIVE);

  const int size = sleqp_working_set_size(working_set);

  working_set->constraint_indices[index] = working_set->num_active_variables
    + (working_set->num_active_constraints)++;

  working_set->content_indices[size] = working_set->num_variables + index;

  working_set->constraint_states[index] = state;

  return SLEQP_OKAY;
}

int sleqp_working_set_get_constraint_index(SleqpWorkingSet* working_set,
                                          int index)
{
  assert(index < working_set->num_constraints);

  return working_set->constraint_indices[index];
}

int sleqp_working_set_get_variable_index(SleqpWorkingSet* working_set,
                                        int index)
{
  assert(index < working_set->num_variables);

  return working_set->variable_indices[index];
}

int sleqp_working_set_get_content(SleqpWorkingSet* working_set,
                                 int index)
{
  return working_set->content_indices[index];
}

SLEQP_ACTIVE_STATE* sleqp_working_set_variable_states(SleqpWorkingSet* working_set)
{
  return working_set->variable_states;
}

SLEQP_ACTIVE_STATE* sleqp_working_set_constraint_states(SleqpWorkingSet* working_set)
{
  return working_set->constraint_states;
}

SLEQP_ACTIVE_STATE sleqp_working_set_get_variable_state(SleqpWorkingSet* working_set,
                                                        int j)
{
  assert(j < working_set->num_variables);
  return working_set->variable_states[j];
}

SLEQP_ACTIVE_STATE sleqp_working_set_get_constraint_state(SleqpWorkingSet* working_set,
                                                          int i)
{
  assert(i < working_set->num_constraints);
  return working_set->constraint_states[i];
}

int sleqp_working_set_num_active_vars(SleqpWorkingSet* working_set)
{
  return working_set->num_active_variables;
}

int sleqp_working_set_num_active_cons(SleqpWorkingSet* working_set)
{
  return working_set->num_active_constraints;
}

int sleqp_working_set_size(SleqpWorkingSet* working_set)
{
  return working_set->num_active_constraints + working_set->num_active_variables;
}

SLEQP_RETCODE sleqp_working_set_fprintf(SleqpWorkingSet* working_set,
                                       FILE* output)
{
  SleqpProblem* problem = working_set->problem;

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
    SLEQP_ACTIVE_STATE state = sleqp_working_set_get_variable_state(working_set, j);

    fprintf(output,
            "State of variable %d: %s\n",
            j,
            state_names[state]);
  }

  for(int i = 0; i < num_constraints; ++i)
  {
    SLEQP_ACTIVE_STATE state = sleqp_working_set_get_constraint_state(working_set, i);

    fprintf(output,
            "State of constraint %d: %s\n",
            i,
            state_names[state]);
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_working_set_copy(SleqpWorkingSet* source,
                                    SleqpWorkingSet* target)
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

SLEQP_RETCODE sleqp_working_set_capture(SleqpWorkingSet* working_set)
{
  ++working_set->refcount;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE working_set_free(SleqpWorkingSet** star)
{
  SleqpWorkingSet* working_set = *star;

  if(!working_set)
  {
    return SLEQP_OKAY;
  }

  sleqp_free(&working_set->content_indices);

  sleqp_free(&working_set->constraint_indices);
  sleqp_free(&working_set->variable_indices);

  sleqp_free(&working_set->constraint_states);
  sleqp_free(&working_set->variable_states);

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_working_set_release(SleqpWorkingSet** star)
{
  SleqpWorkingSet* working_set = *star;

  if(!working_set)
  {
    return SLEQP_OKAY;
  }

  if(--working_set->refcount == 0)
  {
    SLEQP_CALL(working_set_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
