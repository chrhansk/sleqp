#include "sleqp_preprocessing_state.h"

#include "sleqp_mem.h"

struct SleqpPreprocessingState
{
  int refcount;
  SleqpProblem* original_problem;

  SleqpConvertedBound* converted_bounds;
  int num_converted_bounds;

  int num_redundant_constraints;

  int num_fixed_variables;

  SleqpVariableState* var_states;
  SleqpConstraintState* cons_states;
};


SLEQP_RETCODE sleqp_preprocessing_state_create(SleqpPreprocessingState** star,
                                               SleqpProblem* problem)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpPreprocessingState* state = *star;

  *state = (SleqpPreprocessingState) {0};
  state->refcount = 1;

  state->original_problem = problem;
  SLEQP_CALL(sleqp_problem_capture(state->original_problem));

  const int num_variables = sleqp_problem_num_variables(problem);
  const int num_linear = sleqp_problem_num_linear_constraints(problem);

  SLEQP_CALL(sleqp_alloc_array(&state->var_states, num_variables));
  SLEQP_CALL(sleqp_alloc_array(&state->cons_states, num_linear));

  SLEQP_CALL(sleqp_alloc_array(&state->converted_bounds, num_linear));

  SLEQP_CALL(sleqp_preprocessing_state_reset(state));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_preprocessing_state_reset(SleqpPreprocessingState* state)
{
  SleqpProblem* problem = state->original_problem;

  const int num_variables = sleqp_problem_num_variables(problem);
  const int num_linear = sleqp_problem_num_linear_constraints(problem);

  for(int j = 0; j < num_variables; ++j)
  {
    state->var_states[j] = (SleqpVariableState) { .state = SLEQP_VAR_UNCHANGED };
  }

  for(int j = 0; j < num_linear; ++j)
  {
    state->cons_states[j] = (SleqpConstraintState) { .state = SLEQP_CONS_UNCHANGED };
  }

  state->num_converted_bounds = 0;

  state->num_redundant_constraints = 0;

  state->num_fixed_variables = 0;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_preprocessing_state_convert_linear_constraint_to_bound(SleqpPreprocessingState* state,
                                                                           int constraint,
                                                                           int variable,
                                                                           double factor,
                                                                           double var_lb,
                                                                           double var_ub,
                                                                           SleqpBoundState bound_state)
{
  SleqpProblem* problem = state->original_problem;

  const int num_variables = sleqp_problem_num_variables(problem);
  const int num_linear = sleqp_problem_num_linear_constraints(problem);

  assert(constraint >= 0);
  assert(constraint < num_linear);

  assert(variable >= 0);
  assert(variable < num_variables);

  assert(state->var_states[variable].state == SLEQP_VAR_UNCHANGED);
  assert(state->cons_states[constraint].state == SLEQP_CONS_UNCHANGED);

  state->converted_bounds[state->num_converted_bounds] = (SleqpConvertedBound)
    { .constraint = constraint,
      .variable = variable,
      .factor = factor,
      .var_lb = var_lb,
      .var_ub = var_ub,
      .state = bound_state
    };

  state->cons_states[constraint] = (SleqpConstraintState)
    { .state = SLEQP_CONS_BOUNDCONVERTED,
      .bound = state->num_converted_bounds
    };

  ++(state->num_converted_bounds);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_preprocessing_state_remove_linear_constraint(SleqpPreprocessingState* state,
                                                                 int constraint)
{
  SleqpProblem* problem = state->original_problem;

  const int num_linear = sleqp_problem_num_linear_constraints(problem);

  assert(constraint >= 0);
  assert(constraint < num_linear);

  assert(state->cons_states[constraint].state == SLEQP_CONS_UNCHANGED);

  state->cons_states[constraint].state = SLEQP_CONS_REDUNDANT;

  ++(state->num_redundant_constraints);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_preprocessing_state_fix_variable(SleqpPreprocessingState* state,
                                                     int variable,
                                                     double value)
{
  SleqpProblem* problem = state->original_problem;

  const int num_variables = sleqp_problem_num_variables(problem);

  assert(variable >= 0);
  assert(variable < num_variables);

  assert(state->var_states[variable].state == SLEQP_VAR_UNCHANGED);

  state->var_states[variable] = (SleqpVariableState)
  {
    .state = SLEQP_VAR_BOUNDFIXED,
    .value = value
  };

  ++(state->num_fixed_variables);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_preprocessing_state_converted_bounds(SleqpPreprocessingState* state,
                                                         SleqpConvertedBound** star,
                                                         int* num_converted_bounds)
{
  (*star) = state->converted_bounds;

  (*num_converted_bounds) = state->num_converted_bounds;

  return SLEQP_OKAY;
}

SleqpVariableState* sleqp_preprocessing_state_variable_states(const SleqpPreprocessingState* state)
{
  return state->var_states;
}

SleqpConstraintState* sleqp_preprocessing_state_linear_constraint_states(const SleqpPreprocessingState* state)
{
  return state->cons_states;
}

int sleqp_preprocessing_state_num_fixed_variables(const SleqpPreprocessingState* state)
{
  return state->num_fixed_variables;
}

int sleqp_preprocessing_state_num_removed_linear_constraints(const SleqpPreprocessingState* state)
{
  return state->num_redundant_constraints + state->num_converted_bounds;
}

SleqpProblem* sleqp_preprocessing_state_get_problem(const SleqpPreprocessingState* state)
{
  return state->original_problem;
}

static
SLEQP_RETCODE preprocessing_state_free(SleqpPreprocessingState** star)
{
  SleqpPreprocessingState* state = *star;

  sleqp_free(&state->converted_bounds);
  sleqp_free(&state->cons_states);
  sleqp_free(&state->var_states);

  SLEQP_CALL(sleqp_problem_release(&state->original_problem));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_preprocessing_state_capture(SleqpPreprocessingState* state)
{
  ++state->refcount;
  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_preprocessing_state_release(SleqpPreprocessingState** star)
{
  SleqpPreprocessingState* state = *star;

  if(!state)
  {
    return SLEQP_OKAY;
  }

  if(--state->refcount == 0)
  {
    SLEQP_CALL(preprocessing_state_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
