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

  // a map of 0..set_size - 1 -> (variable index) or (num_problem_variables + constraint index)
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
  SleqpProblem* problem = working_set->problem;

  assert(index >= 0);
  assert(index < problem->num_variables);
  assert(working_set->num_active_constraints == 0);

  assert(sleqp_working_set_get_variable_state(working_set, index) == SLEQP_INACTIVE);
  assert(state != SLEQP_INACTIVE);

  if(sleqp_working_set_num_active_cons(working_set) != 0)
  {
    sleqp_log_error("Must add variables before constraints");
    return SLEQP_INTERNAL_ERROR;
  }

  const int size = sleqp_working_set_size(working_set);

  working_set->variable_indices[index] = (working_set->num_active_variables);

  ++(working_set->num_active_variables);

  working_set->content_indices[size] = index;

  working_set->variable_states[index] = state;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_working_set_add_constraint(SleqpWorkingSet* working_set,
                                               int index,
                                               SLEQP_ACTIVE_STATE state)
{
  SleqpProblem* problem = working_set->problem;

  assert(index >= 0);
  assert(index < problem->num_constraints);

  assert(sleqp_working_set_get_constraint_state(working_set, index) == SLEQP_INACTIVE);
  assert(state != SLEQP_INACTIVE);

  const int size = sleqp_working_set_size(working_set);

  working_set->constraint_indices[index] = working_set->num_active_variables
    + (working_set->num_active_constraints);

  ++(working_set->num_active_constraints);

  working_set->content_indices[size] = working_set->num_variables + index;

  working_set->constraint_states[index] = state;

  return SLEQP_OKAY;
}

int sleqp_working_set_get_constraint_index(const SleqpWorkingSet* working_set,
                                           int index)
{
  assert(index < working_set->num_constraints);

  return working_set->constraint_indices[index];
}

int sleqp_working_set_get_variable_index(const SleqpWorkingSet* working_set,
                                         int index)
{
  assert(index < working_set->num_variables);

  return working_set->variable_indices[index];
}

int sleqp_working_set_get_content(const SleqpWorkingSet* working_set,
                                  int index)
{
  return working_set->content_indices[index];
}

SLEQP_ACTIVE_STATE* sleqp_working_set_variable_states(const SleqpWorkingSet* working_set)
{
  return working_set->variable_states;
}

SLEQP_ACTIVE_STATE* sleqp_working_set_constraint_states(const SleqpWorkingSet* working_set)
{
  return working_set->constraint_states;
}

SLEQP_ACTIVE_STATE sleqp_working_set_get_variable_state(const SleqpWorkingSet* working_set,
                                                        int j)
{
  assert(j < working_set->num_variables);
  return working_set->variable_states[j];
}

SLEQP_ACTIVE_STATE sleqp_working_set_get_constraint_state(const SleqpWorkingSet* working_set,
                                                          int i)
{
  assert(i < working_set->num_constraints);
  return working_set->constraint_states[i];
}

int sleqp_working_set_num_active_vars(const SleqpWorkingSet* working_set)
{
  return working_set->num_active_variables;
}

int sleqp_working_set_num_active_cons(const SleqpWorkingSet* working_set)
{
  return working_set->num_active_constraints;
}

int sleqp_working_set_size(const SleqpWorkingSet* working_set)
{
  return working_set->num_active_constraints + working_set->num_active_variables;
}

bool sleqp_working_set_valid(const SleqpWorkingSet* working_set)
{
  SleqpProblem* problem = working_set->problem;

  const int num_variables = problem->num_variables;
  const int num_constraints = problem->num_constraints;

  const int working_set_size = sleqp_working_set_size(working_set);

  {
    int num_active_vars = 0, num_active_cons = 0;

    for(int j = 0; j < num_variables; ++j)
    {
      if(sleqp_working_set_get_variable_state(working_set, j) != SLEQP_INACTIVE)
      {
        ++num_active_vars;
      }
    }

    for(int i = 0; i < num_constraints; ++i)
    {
      if(sleqp_working_set_get_constraint_state(working_set, i) != SLEQP_INACTIVE)
      {
        ++num_active_cons;
      }
    }

    if(num_active_vars != sleqp_working_set_num_active_vars(working_set) ||
       num_active_cons != sleqp_working_set_num_active_cons(working_set))
    {
      return false;
    }
  }

  {
    int num_active_vars = 0;

    for(int j = 0; j < num_variables; ++j)
    {
      int j_idx = sleqp_working_set_get_variable_index(working_set, j);

      if(j_idx == -1)
      {
        continue;
      }

      // variables must appear before constraints
      if(j_idx >= sleqp_working_set_num_active_vars(working_set))
      {
        return false;
      }

      ++num_active_vars;

      for(int k = 0; k < num_variables; ++k)
      {
        int k_idx = sleqp_working_set_get_variable_index(working_set, k);

        // ensure indices are unique
        if((j == k) != (j_idx == k_idx))
        {
          return false;
        }
      }
    }

    if(num_active_vars != sleqp_working_set_num_active_vars(working_set))
    {
      return false;
    }
  }

  {
    int num_active_cons = 0;

    for(int i = 0; i < num_constraints; ++i)
    {
      int i_idx = sleqp_working_set_get_constraint_index(working_set, i);

      if(i_idx == -1)
      {
        continue;
      }

      if(i_idx >= working_set_size)
      {
        return false;
      }

      ++num_active_cons;

      for(int k = 0; k < num_constraints; ++k)
      {
        int k_idx = sleqp_working_set_get_constraint_index(working_set, k);

        // ensure indices are unique
        if((i == k) != (i_idx == k_idx))
        {
          return false;
        }
      }
    }

    if(num_active_cons != sleqp_working_set_num_active_cons(working_set))
    {
      return false;
    }
  }

  for(int k = 0; k < working_set_size; ++k)
  {
    int k_idx = sleqp_working_set_get_content(working_set, k);

    if(k_idx < 0 || k_idx >= (num_variables + num_constraints))
    {
      return false;
    }

    if(k_idx < num_variables)
    {
      const int j = k_idx;

      if(sleqp_working_set_get_variable_state(working_set, j) == SLEQP_INACTIVE)
      {
        return false;
      }
    }
    else
    {
      const int i = k_idx - num_variables;

      if(sleqp_working_set_get_constraint_state(working_set, i) == SLEQP_INACTIVE)
      {
        return false;
      }
    }

    for(int l = 0; l < working_set_size; ++l)
    {
      int l_idx = sleqp_working_set_get_content(working_set, l);

      // ensure indices are unique
      if((l == k) != (l_idx == k_idx))
      {
        return false;
      }
    }
  }

  return true;
}

SLEQP_RETCODE sleqp_working_set_fprintf(const SleqpWorkingSet* working_set,
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

SLEQP_RETCODE sleqp_working_set_copy(const SleqpWorkingSet* source,
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
