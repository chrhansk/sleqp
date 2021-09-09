#include "working_set.h"

#include <assert.h>

#include "fail.h"
#include "mem.h"

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

  // a mapping of 0..num_variables - 1 -> pos in the active set or SLEQP_NONE
  int* variable_indices;

  // a mapping of 0..num_constraints - 1 -> pos in the active set or SLEQP_NONE
  int* constraint_indices;

  // a map of 0..set_size - 1 -> (variable index) or (num_problem_variables + constraint index)
  int* content_indices;

};

SLEQP_RETCODE sleqp_working_set_create(SleqpWorkingSet** star,
                                       SleqpProblem* problem)
{
  SLEQP_CALL(sleqp_malloc(star));

  const int num_variables = sleqp_problem_num_variables(problem);
  const int num_constraints = sleqp_problem_num_constraints(problem);

  SleqpWorkingSet* working_set = *star;

  working_set->refcount = 1;

  working_set->problem = problem;
  SLEQP_CALL(sleqp_problem_capture(working_set->problem));

  SLEQP_CALL(sleqp_alloc_array(&working_set->variable_states, num_variables));
  SLEQP_CALL(sleqp_alloc_array(&working_set->constraint_states, num_constraints));

  working_set->num_variables = num_variables;
  working_set->num_constraints = num_constraints;

  working_set->max_set_size = num_variables + num_constraints;

  SLEQP_CALL(sleqp_alloc_array(&working_set->variable_indices, num_variables));
  SLEQP_CALL(sleqp_alloc_array(&working_set->constraint_indices, num_constraints));

  SLEQP_CALL(sleqp_alloc_array(&working_set->content_indices,
                               working_set->max_set_size));

  SLEQP_CALL(sleqp_working_set_reset(working_set));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_working_set_reset(SleqpWorkingSet* working_set)
{
  for(int j = 0; j < working_set->num_variables; ++j)
  {
    working_set->variable_states[j] = SLEQP_INACTIVE;
    working_set->variable_indices[j] = SLEQP_NONE;
  }

  for(int i = 0; i < working_set->num_constraints; ++i)
  {
    working_set->constraint_states[i] = SLEQP_INACTIVE;
    working_set->constraint_indices[i] = SLEQP_NONE;
  }

  working_set->num_active_constraints = 0;
  working_set->num_active_variables = 0;

  return SLEQP_OKAY;
}

bool sleqp_working_set_eq(SleqpWorkingSet* first,
                          SleqpWorkingSet* second)
{
  assert(first->num_constraints == second->num_constraints);
  assert(first->num_variables == second->num_variables);

  for(int j = 0; j < first->num_variables; ++j)
  {
    if(first->variable_indices[j] != second->variable_indices[j])
    {
      return false;
    }
  }

  for(int i = 0; i < first->num_constraints; ++i)
  {
    if(first->constraint_indices[i] != second->constraint_indices[i])
    {
      return false;
    }
  }

  return true;
}

SLEQP_RETCODE sleqp_working_set_add_variable(SleqpWorkingSet* working_set,
                                             int index,
                                             SLEQP_ACTIVE_STATE state)
{
  SleqpProblem* problem = working_set->problem;

  const int num_variables = sleqp_problem_num_variables(problem);

  assert(index >= 0);
  assert(index < num_variables);
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

  const int num_constraints = sleqp_problem_num_constraints(problem);

  assert(index >= 0);
  assert(index < num_constraints);

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

SLEQP_RETCODE sleqp_working_set_add(SleqpWorkingSet* working_set,
                                    int index,
                                    bool constraint,
                                    SLEQP_ACTIVE_STATE state)
{
  if(constraint)
  {
    return sleqp_working_set_add_constraint(working_set, index, state);
  }
  else
  {
    return sleqp_working_set_add_variable(working_set, index, state);
  }
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

SLEQP_ACTIVE_STATE sleqp_working_set_get_state(const SleqpWorkingSet* working_set,
                                               bool constraint,
                                               int index)
{
  if(constraint)
  {
    return sleqp_working_set_get_constraint_state(working_set, index);
  }
  else
  {
    return sleqp_working_set_get_variable_state(working_set, index);
  }
}

SleqpProblem* sleqp_working_set_get_problem(const SleqpWorkingSet* working_set)
{
  return working_set->problem;
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

  const int num_variables = sleqp_problem_num_variables(problem);
  const int num_constraints = sleqp_problem_num_constraints(problem);

  const int working_set_size = sleqp_working_set_size(working_set);

  {
    const int num_active_vars = sleqp_working_set_num_active_vars(working_set);
    const int num_active_cons = sleqp_working_set_num_active_cons(working_set);
    const int working_set_size = sleqp_working_set_size(working_set);

    assert(num_active_vars <= working_set_size);
    assert(num_active_cons <= working_set_size);

    assert(working_set_size <= num_variables);
  }

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

      if(j_idx == SLEQP_NONE)
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

      if(i_idx == SLEQP_NONE)
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

static
SLEQP_RETCODE supports_range(SleqpSparseVec* dual,
                             SLEQP_ACTIVE_STATE* states,
                             bool* supports)
{
  *supports = true;

  for(int k = 0; k < dual->nnz; ++k)
  {
    const int i = dual->indices[k];
    const double v = dual->data[k];

    SLEQP_ACTIVE_STATE state = states[i];

    switch (state)
    {
    case SLEQP_ACTIVE_BOTH:
      break;
    case SLEQP_ACTIVE_LOWER:
      if(v > 0.)
      {
        *supports = false;
        return SLEQP_OKAY;
      }
      break;
    case SLEQP_ACTIVE_UPPER:
      if(v < 0.)
      {
        *supports = false;
        return SLEQP_OKAY;
      }
      break;
    case SLEQP_INACTIVE:
      if(v != 0.)
      {
        *supports = false;
        return SLEQP_OKAY;
      }
      break;
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_working_set_supports_cons_dual(const SleqpWorkingSet* working_set,
                                                   SleqpSparseVec* cons_dual,
                                                   bool* supports)
{
  SleqpProblem* problem = working_set->problem;

  sleqp_assert(cons_dual->dim == sleqp_problem_num_constraints(problem));

  SLEQP_CALL(supports_range(cons_dual, working_set->constraint_states, supports));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_working_set_supports_vars_dual(const SleqpWorkingSet* working_set,
                                                   SleqpSparseVec* vars_dual,
                                                   bool* supports)
{
  SleqpProblem* problem = working_set->problem;

  sleqp_assert(vars_dual->dim == sleqp_problem_num_variables(problem));

  SLEQP_CALL(supports_range(vars_dual, working_set->variable_states, supports));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_working_set_fprintf(const SleqpWorkingSet* working_set,
                                        FILE* output)
{
  SleqpProblem* problem = working_set->problem;

  const int num_variables = sleqp_problem_num_variables(problem);
  const int num_constraints = sleqp_problem_num_constraints(problem);

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
  SleqpProblem* source_problem = source->problem;
  SleqpProblem* target_problem = target->problem;

  const int num_variables = sleqp_problem_num_variables(source_problem);
  const int num_constraints = sleqp_problem_num_constraints(source_problem);

  const int max_set_size = source->max_set_size;

  assert(num_variables == sleqp_problem_num_variables(target_problem));
  assert(num_constraints == sleqp_problem_num_constraints(target_problem));
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

  SLEQP_CALL(sleqp_problem_release(&working_set->problem));

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
