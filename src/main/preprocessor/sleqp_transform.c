#include "sleqp_transform.h"

#include "sleqp_assert.h"
#include "sleqp_cmp.h"
#include "sleqp_mem.h"

#include "preprocessor/sleqp_preprocessing.h"

#include "preprocessor/sleqp_fixed_var_func.h"

struct SleqpTransformation
{
  int refcount;

  int num_transformed_vars;
  int num_transformed_linear_cons;

  SleqpPreprocessingState* preprocessing_state;
  SleqpProblem* original_problem;

  SleqpParams* params;

  int num_fixed_vars;
  int* fixed_var_indices;
  double* fixed_var_values;

  int num_removed_cons;
  int* removed_cons_indices;

  SleqpSparseVec* transformed_var_lb;
  SleqpSparseVec* transformed_var_ub;

  SleqpSparseVec* transformed_linear_lb;
  SleqpSparseVec* transformed_linear_ub;

  SleqpSparseMatrix* transformed_linear_coeffs;

  SleqpFunc* transformed_func;

  double* dense_cache;
  SleqpSparseVec* sparse_cache;
};

static
SLEQP_RETCODE create_fixed_variables(SleqpTransformation* transformation)
{
  SleqpPreprocessingState* preprocessing_state = transformation->preprocessing_state;

  SleqpProblem* problem = transformation->original_problem;

  const int num_variables = sleqp_problem_num_variables(problem);

  transformation->num_fixed_vars = sleqp_preprocessing_state_num_fixed_variables(preprocessing_state);

  SLEQP_CALL(sleqp_alloc_array(&transformation->fixed_var_indices,
                               transformation->num_fixed_vars));

  SLEQP_CALL(sleqp_alloc_array(&transformation->fixed_var_values,
                               transformation->num_fixed_vars));

  SleqpVariableState* var_states = sleqp_preprocessing_state_variable_states(preprocessing_state);

  int i = 0;

  for(int j = 0; j < num_variables; ++j)
  {
    if(var_states[j].state != SLEQP_VAR_UNCHANGED)
    {
      assert(var_states[j].state == SLEQP_VAR_BOUNDFIXED);

      transformation->fixed_var_indices[i] = j;
      transformation->fixed_var_values[i] = var_states[j].value;

      ++i;
    }
  }

  assert(i == transformation->num_fixed_vars);

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE create_fixed_linear_constraints(SleqpTransformation* transformation)
{
  SleqpPreprocessingState* preprocessing_state = transformation->preprocessing_state;

  SleqpProblem* problem = transformation->original_problem;

  const int num_linear_cons = sleqp_problem_num_linear_constraints(problem);

   transformation->num_removed_cons = sleqp_preprocessing_state_num_removed_linear_constraints(preprocessing_state);

  SLEQP_CALL(sleqp_alloc_array(&transformation->removed_cons_indices, transformation->num_removed_cons));

  SleqpConstraintState* cons_states = sleqp_preprocessing_state_linear_constraint_states(preprocessing_state);

  int j = 0;

  for(int i = 0; i < num_linear_cons; ++i)
  {
    if(cons_states[i].state != SLEQP_CONS_UNCHANGED)
    {
      transformation->removed_cons_indices[j] = i;
      ++j;
    }
  }

  assert(j == transformation->num_removed_cons);

  return SLEQP_OKAY;
}


SLEQP_RETCODE sleqp_transformation_create(SleqpTransformation** star,
                                          SleqpPreprocessingState* preprocessing_state,
                                          SleqpParams* params)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpTransformation* transformation = *star;

  *transformation = (SleqpTransformation) {0};
  transformation->refcount = 1;

  SleqpProblem* problem = sleqp_preprocessing_state_get_problem(preprocessing_state);

  transformation->original_problem = problem;
  SLEQP_CALL(sleqp_problem_capture(transformation->original_problem));

  transformation->params = params;
  SLEQP_CALL(sleqp_params_capture(transformation->params));

  transformation->preprocessing_state = preprocessing_state;
  SLEQP_CALL(sleqp_preprocessing_state_capture(transformation->preprocessing_state));

  SLEQP_CALL(create_fixed_variables(transformation));

  SLEQP_CALL(create_fixed_linear_constraints(transformation));

  const int num_variables = sleqp_problem_num_variables(problem);

  assert(transformation->num_fixed_vars >= 0);
  assert(transformation->num_fixed_vars < num_variables);

  transformation->num_transformed_vars = num_variables - transformation->num_fixed_vars;

  const int num_linear_cons = sleqp_problem_num_linear_constraints(problem);

  const int num_removed_cons = sleqp_preprocessing_state_num_removed_linear_constraints(preprocessing_state);

  transformation->num_transformed_linear_cons = num_linear_cons - num_removed_cons;

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&transformation->transformed_var_lb,
                                              transformation->num_transformed_vars));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&transformation->transformed_var_ub,
                                              transformation->num_transformed_vars));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&transformation->transformed_linear_lb,
                                              transformation->num_transformed_linear_cons));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&transformation->transformed_linear_ub,
                                              transformation->num_transformed_linear_cons));

  SLEQP_CALL(sleqp_sparse_matrix_create(&transformation->transformed_linear_coeffs,
                                        transformation->num_transformed_linear_cons,
                                        transformation->num_transformed_vars,
                                        0));

  SLEQP_CALL(sleqp_alloc_array(&transformation->dense_cache, num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_full(&transformation->sparse_cache,
                                             transformation->num_transformed_vars));

  return SLEQP_OKAY;
}


SLEQP_RETCODE sleqp_transformation_convert_primal(SleqpTransformation* transformation,
                                                  const SleqpSparseVec* source,
                                                  SleqpSparseVec* target)
{
  SLEQP_CALL(sleqp_preprocessing_remove_entries(source,
                                                target,
                                                transformation->num_fixed_vars,
                                                transformation->fixed_var_indices));

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE create_transformed_func(SleqpTransformation* transformation,
                                      SleqpFunc** star)
{
  SleqpProblem* problem = transformation->original_problem;
  SleqpFunc* func = sleqp_problem_func(problem);

  if(transformation->num_fixed_vars == 0)
  {
    SLEQP_CALL(sleqp_func_capture(func));

    *star = func;
  }
  else
  {
    SLEQP_CALL(sleqp_fixed_var_func_create(star,
                                           func,
                                           transformation->num_fixed_vars,
                                           transformation->fixed_var_indices,
                                           transformation->fixed_var_values));
  }

  assert(sleqp_func_get_num_variables(*star) ==
         sleqp_problem_num_variables(problem) - transformation->num_fixed_vars);

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE create_transformed_var_lb(SleqpTransformation* transformation)
{
  SleqpProblem* problem = transformation->original_problem;

  const int num_variables = sleqp_problem_num_variables(problem);

  const double zero_eps = sleqp_params_get(transformation->params,
                                           SLEQP_PARAM_ZERO_EPS);

  SLEQP_CALL(sleqp_sparse_vector_to_raw(sleqp_problem_var_lb(problem),
                                        transformation->dense_cache));

  SleqpConvertedBound* converted_bounds;
  int num_converted_bounds;

  SLEQP_CALL(sleqp_preprocessing_state_converted_bounds(transformation->preprocessing_state,
                                                        &converted_bounds,
                                                        &num_converted_bounds));

  for(int k = 0; k < num_converted_bounds; ++k)
  {
    const int i = converted_bounds[k].constraint;

    transformation->dense_cache[i] = SLEQP_MAX(transformation->dense_cache[i],
                                               converted_bounds[k].var_lb);
  }

  SLEQP_CALL(sleqp_sparse_vector_from_raw(transformation->sparse_cache,
                                          transformation->dense_cache,
                                          num_variables,
                                          zero_eps));

  SLEQP_CALL(sleqp_preprocessing_remove_entries(transformation->sparse_cache,
                                                transformation->transformed_var_lb,
                                                transformation->num_fixed_vars,
                                                transformation->fixed_var_indices));

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE create_transformed_var_ub(SleqpTransformation* transformation)
{
  SleqpProblem* problem = transformation->original_problem;

  const int num_variables = sleqp_problem_num_variables(problem);

  const double zero_eps = sleqp_params_get(transformation->params,
                                           SLEQP_PARAM_ZERO_EPS);

  SLEQP_CALL(sleqp_sparse_vector_to_raw(sleqp_problem_var_ub(problem),
                                        transformation->dense_cache));

  SleqpConvertedBound* converted_bounds;
  int num_converted_bounds;

  SLEQP_CALL(sleqp_preprocessing_state_converted_bounds(transformation->preprocessing_state,
                                                        &converted_bounds,
                                                        &num_converted_bounds));

  for(int k = 0; k < num_converted_bounds; ++k)
  {
    const int i = converted_bounds[k].constraint;

    transformation->dense_cache[i] = SLEQP_MIN(transformation->dense_cache[i],
                                               converted_bounds[k].var_ub);
  }

  SLEQP_CALL(sleqp_sparse_vector_from_raw(transformation->sparse_cache,
                                          transformation->dense_cache,
                                          num_variables,
                                          zero_eps));

  SLEQP_CALL(sleqp_preprocessing_remove_entries(transformation->sparse_cache,
                                                transformation->transformed_var_ub,
                                                transformation->num_fixed_vars,
                                                transformation->fixed_var_indices));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_transformation_create_transformed_problem(SleqpTransformation* transformation,
                                                              SleqpProblem** star)
{
  SleqpProblem* problem = transformation->original_problem;

  SLEQP_CALL(create_transformed_func(transformation, &transformation->transformed_func));

  SLEQP_CALL(create_transformed_var_lb(transformation));

  SLEQP_CALL(create_transformed_var_ub(transformation));

  SLEQP_CALL(sleqp_preprocessing_remove_entries(sleqp_problem_linear_lb(problem),
                                                transformation->transformed_linear_lb,
                                                transformation->num_removed_cons,
                                                transformation->removed_cons_indices));

  SLEQP_CALL(sleqp_preprocessing_remove_entries(sleqp_problem_linear_ub(problem),
                                                transformation->transformed_linear_ub,
                                                transformation->num_removed_cons,
                                                transformation->removed_cons_indices));

  SLEQP_CALL(sleqp_preprocessing_remove_matrix_entries(sleqp_problem_linear_coeffs(problem),
                                                       transformation->transformed_linear_coeffs,
                                                       transformation->num_fixed_vars,
                                                       transformation->fixed_var_indices,
                                                       transformation->num_removed_cons,
                                                       transformation->removed_cons_indices));

  SLEQP_CALL(sleqp_problem_create(star,
                                  transformation->transformed_func,
                                  transformation->params,
                                  transformation->transformed_var_lb,
                                  transformation->transformed_var_ub,
                                  sleqp_problem_general_lb(problem),
                                  sleqp_problem_general_ub(problem),
                                  transformation->transformed_linear_coeffs,
                                  transformation->transformed_linear_lb,
                                  transformation->transformed_linear_ub));

  {
    SleqpProblem* transformed_problem = *star;

    const int num_variables = sleqp_problem_num_variables(problem);
    const int num_constraints = sleqp_problem_num_constraints(problem);

    const int num_transformed_variables = sleqp_problem_num_variables(transformed_problem);
    const int num_transformed_constraints = sleqp_problem_num_constraints(transformed_problem);

    sleqp_assert(num_transformed_variables + transformation->num_fixed_vars == num_variables);
    sleqp_assert(num_transformed_constraints + transformation->num_removed_cons == num_constraints);
  }

  SLEQP_CALL(sleqp_func_release(&transformation->transformed_func));

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE transformation_free(SleqpTransformation** star)
{
  SleqpTransformation* transformation = *star;

  SLEQP_CALL(sleqp_sparse_vector_free(&transformation->sparse_cache));

  sleqp_free(&transformation->dense_cache);

  SLEQP_CALL(sleqp_func_release(&transformation->transformed_func));

  SLEQP_CALL(sleqp_sparse_matrix_release(&transformation->transformed_linear_coeffs));

  SLEQP_CALL(sleqp_sparse_vector_free(&transformation->transformed_linear_ub));
  SLEQP_CALL(sleqp_sparse_vector_free(&transformation->transformed_linear_lb));

  SLEQP_CALL(sleqp_sparse_vector_free(&transformation->transformed_var_ub));
  SLEQP_CALL(sleqp_sparse_vector_free(&transformation->transformed_var_lb));

  sleqp_free(&transformation->fixed_var_values);
  sleqp_free(&transformation->fixed_var_indices);

  sleqp_free(&transformation->removed_cons_indices);

  SLEQP_CALL(sleqp_preprocessing_state_release(&transformation->preprocessing_state));
  SLEQP_CALL(sleqp_params_release(&transformation->params));
  SLEQP_CALL(sleqp_problem_release(&transformation->original_problem));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_transformation_capture(SleqpTransformation* state)
{
  ++state->refcount;
  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_transformation_release(SleqpTransformation** star)
{
  SleqpTransformation* transformation = *star;

  if(!transformation)
  {
    return SLEQP_OKAY;
  }

  if(--transformation->refcount == 0)
  {
    SLEQP_CALL(transformation_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
