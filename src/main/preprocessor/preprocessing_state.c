#include "preprocessing_state.h"

#include "cmp.h"
#include "log.h"
#include "mem.h"

struct SleqpPreprocessingState
{
  int refcount;
  SleqpProblem* original_problem;

  SleqpConvertedBound* converted_bounds;
  int num_converted_bounds;

  SleqpForcingConstraint* forcing_constraints;
  int num_forcing_constraints;

  int num_redundant_constraints;

  SleqpVariableState* var_states;
  SleqpConstraintState* cons_states;

  SleqpBoundRequirementState* var_bound_states;
  SleqpBoundRequirementState* cons_bound_states;

  int num_fixed_vars;
  int* fixed_var_indices;
  double* fixed_var_values;

  int num_removed_cons;
  int* removed_cons_indices;

  int num_removed_bounds;
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

  SLEQP_CALL(sleqp_alloc_array(&state->var_bound_states, num_variables));
  SLEQP_CALL(sleqp_alloc_array(&state->cons_bound_states, num_linear));

  SLEQP_CALL(sleqp_alloc_array(&state->converted_bounds, num_linear));

  SLEQP_CALL(sleqp_alloc_array(&state->forcing_constraints, num_linear));

  for(int i = 0; i < num_linear; ++i)
  {
    state->forcing_constraints[i] = (SleqpForcingConstraint) {0};
  }

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
    state->var_bound_states[j] = SLEQP_BOUND_REQUIRED;
  }

  for(int j = 0; j < num_linear; ++j)
  {
    state->cons_states[j] = (SleqpConstraintState) { .state = SLEQP_CONS_UNCHANGED };
    state->cons_bound_states[j] = SLEQP_BOUND_REQUIRED;
  }

  state->num_converted_bounds = 0;

  state->num_redundant_constraints = 0;

  state->num_fixed_vars = 0;
  state->num_removed_cons = 0;

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
    {
      .constraint = constraint,
      .variable = variable,
      .factor = factor,
      .var_lb = var_lb,
      .var_ub = var_ub,
      .state = bound_state
    };

  state->cons_states[constraint] = (SleqpConstraintState)
    {
      .state = SLEQP_CONS_BOUNDCONVERTED,
      .bound = state->num_converted_bounds
    };

  ++(state->num_converted_bounds);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_preprocessing_state_add_forcing_constraint(SleqpPreprocessingState* state,
                                                               int constraint,
                                                               SleqpBoundState bound_state,
                                                               double* var_lb,
                                                               double* var_ub)
{
  SleqpProblem* problem = state->original_problem;

  const int num_variables = sleqp_problem_num_variables(problem);
  const int num_linear = sleqp_problem_num_linear_constraints(problem);

  assert(constraint >= 0);
  assert(constraint < num_linear);

  assert(bound_state != SLEQP_BOTH_BOUNDS);

  SleqpVariableState* var_states = sleqp_preprocessing_state_variable_states(state);
  SleqpConstraintState* linear_states = sleqp_preprocessing_state_linear_constraint_states(state);

  assert(linear_states[constraint].state == SLEQP_CONS_UNCHANGED);

  SleqpSparseMatrix* linear_coeffs = sleqp_problem_linear_coeffs(problem);

  double* linear_data = sleqp_sparse_matrix_get_data(linear_coeffs);
  int* linear_rows = sleqp_sparse_matrix_get_rows(linear_coeffs);
  int* linear_cols = sleqp_sparse_matrix_get_cols(linear_coeffs);

  int num_coeffs = 0;

  for(int col = 0; col < num_variables; ++col)
  {
    if(var_states[col].state != SLEQP_VAR_UNCHANGED)
    {
      continue;
    }

    for (int k = linear_cols[col]; k < linear_cols[col + 1]; ++k)
    {
      const int row = linear_rows[k];
      const double value = linear_data[k];

      if(row < constraint)
      {
        continue;
      }
      else if(row > constraint)
      {
        break;
      }

      if(value == 0.)
      {
        continue;
      }

      ++num_coeffs;
    }
  }

  SLEQP_CALL(sleqp_preprocessing_state_remove_linear_constraint(state,
                                                                constraint));

  if(num_coeffs == 0)
  {
    return SLEQP_OKAY;
  }

  SleqpForcingConstraint* forcing_constraint = state->forcing_constraints + state->num_forcing_constraints;

  SLEQP_CALL(sleqp_alloc_array(&forcing_constraint->variables, num_coeffs));
  SLEQP_CALL(sleqp_alloc_array(&forcing_constraint->factors, num_coeffs));

  num_coeffs = 0;

  for(int col = 0; col < num_variables; ++col)
  {
    if(var_states[col].state != SLEQP_VAR_UNCHANGED)
    {
      continue;
    }

    for (int k = linear_cols[col]; k < linear_cols[col + 1]; ++k)
    {
      const int row = linear_rows[k];
      const double coeff_value = linear_data[k];

      if(row < constraint)
      {
        continue;
      }
      else if(row > constraint)
      {
        break;
      }

      if(coeff_value == 0.)
      {
        break;
      }

      forcing_constraint->variables[num_coeffs] = col;
      forcing_constraint->factors[num_coeffs] = coeff_value;

      double fixed_value;

      bool pos_coeff = (coeff_value > 0);
      bool fixed_lower = (bound_state == SLEQP_LOWER_BOUND);

      if(pos_coeff != fixed_lower)
      {
        fixed_value = var_lb[col];
      }
      else
      {
        fixed_value = var_ub[col];
      }

      assert(sleqp_is_finite(fixed_value));

      state->var_states[col] = (SleqpVariableState)
        {
          .state = SLEQP_VAR_FORCING_FIXED,
          .value = fixed_value
        };

      ++(state->num_fixed_vars);

      ++num_coeffs;
    }
  }

  forcing_constraint->num_variables = num_coeffs;
  forcing_constraint->constraint = constraint;
  forcing_constraint->state = bound_state;

  (++state->num_forcing_constraints);

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
  state->cons_bound_states[constraint] = SLEQP_BOUND_REDUNDANT;

  ++(state->num_redundant_constraints);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_preprocessing_state_fix_variable_to_bounds(SleqpPreprocessingState* state,
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
    .state = SLEQP_VAR_BOUND_FIXED,
    .value = value
  };

  ++(state->num_fixed_vars);

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

SLEQP_RETCODE sleqp_preprocessing_state_forcing_constraints(SleqpPreprocessingState* state,
                                                            SleqpForcingConstraint** star,
                                                            int* num_forcing_constraints)
{
  (*star) = state->forcing_constraints;

  (*num_forcing_constraints) = state->num_forcing_constraints;

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

SleqpBoundRequirementState* sleqp_preprocessing_state_variable_bound_requirements(const SleqpPreprocessingState* state)
{
  return state->var_bound_states;
}

SleqpBoundRequirementState* sleqp_preprocessing_state_linear_constraint_bound_requirements(const SleqpPreprocessingState* state)
{
  return state->cons_bound_states;
}

SLEQP_RETCODE sleqp_preprocessing_state_add_variable_bound_requirement(SleqpPreprocessingState* state,
                                                                       int j,
                                                                       SleqpBoundRequirementState requirement_state)
{
  if(state->var_bound_states[j] == SLEQP_BOUND_REQUIRED)
  {
    ++(state->num_removed_bounds);
  }

  state->var_bound_states[j] |= requirement_state;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_preprocessing_state_add_linear_constraint_bound_requirement(SleqpPreprocessingState* state,
                                                                                int i,
                                                                                SleqpBoundRequirementState requirement_state)
{
  if(state->cons_bound_states[i] == SLEQP_BOUND_REQUIRED)
  {
    ++(state->num_removed_bounds);
  }

  state->cons_bound_states[i] |= requirement_state;

  return SLEQP_OKAY;
}

int sleqp_preprocessing_state_num_fixed_variables(const SleqpPreprocessingState* state)
{
  return state->num_fixed_vars;
}

int sleqp_preprocessing_state_num_removed_linear_constraints(const SleqpPreprocessingState* state)
{
  return state->num_redundant_constraints + state->num_converted_bounds;
}

int sleqp_preprocessing_state_num_removed_bounds(const SleqpPreprocessingState* state)
{
  return state->num_removed_bounds;
}

static
SLEQP_RETCODE create_fixed_variables(SleqpPreprocessingState* state)
{
  SleqpProblem* problem = state->original_problem;

  const int num_variables = sleqp_problem_num_variables(problem);

  sleqp_free(&state->fixed_var_indices);
  sleqp_free(&state->fixed_var_values);

  SLEQP_CALL(sleqp_alloc_array(&state->fixed_var_indices,
                               state->num_fixed_vars));

  SLEQP_CALL(sleqp_alloc_array(&state->fixed_var_values,
                               state->num_fixed_vars));

  SleqpVariableState* var_states = sleqp_preprocessing_state_variable_states(state);

  int i = 0;

  for(int j = 0; j < num_variables; ++j)
  {
    if(var_states[j].state != SLEQP_VAR_UNCHANGED)
    {
      assert(var_states[j].state & SLEQP_VAR_FIXED);

      state->fixed_var_indices[i] = j;
      state->fixed_var_values[i] = var_states[j].value;

      ++i;
    }
  }

  assert(i == state->num_fixed_vars);

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE create_removed_linear_constraints(SleqpPreprocessingState* state)
{
  SleqpProblem* problem = state->original_problem;

  const int num_linear_cons = sleqp_problem_num_linear_constraints(problem);

  state->num_removed_cons = sleqp_preprocessing_state_num_removed_linear_constraints(state);

  sleqp_free(&state->removed_cons_indices);

  SLEQP_CALL(sleqp_alloc_array(&state->removed_cons_indices, state->num_removed_cons));

  SleqpConstraintState* cons_states = sleqp_preprocessing_state_linear_constraint_states(state);

  int j = 0;

  for(int i = 0; i < num_linear_cons; ++i)
  {
    if(cons_states[i].state != SLEQP_CONS_UNCHANGED)
    {
      state->removed_cons_indices[j] = i;
      ++j;
    }
  }

  assert(j == state->num_removed_cons);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_preprocessing_state_flush(SleqpPreprocessingState* state)
{
  SLEQP_CALL(create_fixed_variables(state));

  SLEQP_CALL(create_removed_linear_constraints(state));

  SLEQP_CALL(sleqp_realloc(&state->converted_bounds, state->num_converted_bounds));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_preprocessing_state_fixed_variables(SleqpPreprocessingState* state,
                                                        int* num_fixed_vars,
                                                        int** fixed_var_indices,
                                                        double** fixed_var_values)
{
  *num_fixed_vars = state->num_fixed_vars;

  *fixed_var_indices = state->fixed_var_indices;
  *fixed_var_values = state->fixed_var_values;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_preprocessing_state_removed_linear_constraints(SleqpPreprocessingState* state,
                                                                   int* num_removed_cons,
                                                                   int** removed_cons_indices)
{
  *num_removed_cons = sleqp_preprocessing_state_num_removed_linear_constraints(state);

  *removed_cons_indices = state->removed_cons_indices;

  return SLEQP_OKAY;
}

SleqpProblem* sleqp_preprocessing_state_get_problem(const SleqpPreprocessingState* state)
{
  return state->original_problem;
}

static
SLEQP_RETCODE preprocessing_state_free(SleqpPreprocessingState** star)
{
  SleqpPreprocessingState* state = *star;

  sleqp_free(&state->removed_cons_indices);

  sleqp_free(&state->fixed_var_indices);
  sleqp_free(&state->fixed_var_values);

  for(int k = 0; k < state->num_forcing_constraints; ++k)
  {
    sleqp_free(&state->forcing_constraints[k].factors);
    sleqp_free(&state->forcing_constraints[k].variables);
  }

  sleqp_free(&state->forcing_constraints);

  sleqp_free(&state->converted_bounds);

  sleqp_free(&state->cons_bound_states);
  sleqp_free(&state->var_bound_states);

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
