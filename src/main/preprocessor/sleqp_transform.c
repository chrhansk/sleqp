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

  SleqpSparseVec* transformed_var_lb;
  SleqpSparseVec* transformed_var_ub;

  SleqpSparseVec* fixed_variable_values;

  SleqpSparseVec* fixed_linear_lb;
  SleqpSparseVec* fixed_linear_ub;

  SleqpSparseVec* transformed_linear_lb;
  SleqpSparseVec* transformed_linear_ub;

  SleqpSparseMatrix* transformed_linear_coeffs;

  SleqpFunc* transformed_func;

  double* dense_cache;
  SleqpSparseVec* sparse_cache;
};


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

  const int num_variables = sleqp_problem_num_variables(problem);
  const int num_constraints = sleqp_problem_num_constraints(problem);

  const int num_linear = sleqp_problem_num_linear_constraints(problem);

  const int num_fixed_vars = sleqp_preprocessing_state_num_fixed_variables(preprocessing_state);

  transformation->num_transformed_vars = num_variables - num_fixed_vars;

  const int num_linear_cons = sleqp_problem_num_linear_constraints(problem);

  const int num_removed_cons = sleqp_preprocessing_state_num_removed_linear_constraints(preprocessing_state);

  transformation->num_transformed_linear_cons = num_linear_cons - num_removed_cons;

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&transformation->transformed_var_lb,
                                              transformation->num_transformed_vars));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&transformation->transformed_var_ub,
                                              transformation->num_transformed_vars));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&transformation->fixed_variable_values,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&transformation->fixed_linear_lb,
                                              num_linear));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&transformation->fixed_linear_ub,
                                              num_linear));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&transformation->transformed_linear_lb,
                                              transformation->num_transformed_linear_cons));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&transformation->transformed_linear_ub,
                                              transformation->num_transformed_linear_cons));

  SLEQP_CALL(sleqp_sparse_matrix_create(&transformation->transformed_linear_coeffs,
                                        transformation->num_transformed_linear_cons,
                                        transformation->num_transformed_vars,
                                        0));

  SLEQP_CALL(sleqp_alloc_array(&transformation->dense_cache, SLEQP_MAX(num_variables, num_constraints)));

  SLEQP_CALL(sleqp_sparse_vector_create_full(&transformation->sparse_cache,
                                             transformation->num_transformed_vars));

  return SLEQP_OKAY;
}


SLEQP_RETCODE sleqp_transformation_convert_primal(SleqpTransformation* transformation,
                                                  const SleqpSparseVec* source,
                                                  SleqpSparseVec* target)
{
  int num_fixed_vars;
  int* fixed_var_indices;
  double* fixed_var_values;

  SLEQP_CALL(sleqp_preprocessing_state_fixed_variables(transformation->preprocessing_state,
                                                       &num_fixed_vars,
                                                       &fixed_var_indices,
                                                       &fixed_var_values));

  SLEQP_CALL(sleqp_preprocessing_remove_entries(source,
                                                target,
                                                num_fixed_vars,
                                                fixed_var_indices));

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE create_transformed_func(SleqpTransformation* transformation,
                                      SleqpFunc** star)
{
  SleqpProblem* problem = transformation->original_problem;
  SleqpFunc* func = sleqp_problem_func(problem);

  int num_fixed_vars;
  int* fixed_var_indices;
  double* fixed_var_values;

  SLEQP_CALL(sleqp_preprocessing_state_fixed_variables(transformation->preprocessing_state,
                                                       &num_fixed_vars,
                                                       &fixed_var_indices,
                                                       &fixed_var_values));

  if(num_fixed_vars == 0)
  {
    SLEQP_CALL(sleqp_func_capture(func));

    *star = func;
  }
  else
  {
    SLEQP_CALL(sleqp_fixed_var_func_create(star,
                                           func,
                                           num_fixed_vars,
                                           fixed_var_indices,
                                           fixed_var_values));
  }

  assert(sleqp_func_get_num_variables(*star) == sleqp_problem_num_variables(problem) - num_fixed_vars);

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE create_transformed_var_lb(SleqpTransformation* transformation)
{
  SleqpProblem* problem = transformation->original_problem;

  const int num_variables = sleqp_problem_num_variables(problem);

  const double zero_eps = sleqp_params_get(transformation->params,
                                           SLEQP_PARAM_ZERO_EPS);

  int num_fixed_vars;
  int* fixed_var_indices;
  double* fixed_var_values;

  SLEQP_CALL(sleqp_preprocessing_state_fixed_variables(transformation->preprocessing_state,
                                                       &num_fixed_vars,
                                                       &fixed_var_indices,
                                                       &fixed_var_values));

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
                                                num_fixed_vars,
                                                fixed_var_indices));

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE create_transformed_var_ub(SleqpTransformation* transformation)
{
  SleqpProblem* problem = transformation->original_problem;

  const int num_variables = sleqp_problem_num_variables(problem);

  const double zero_eps = sleqp_params_get(transformation->params,
                                           SLEQP_PARAM_ZERO_EPS);

  SleqpPreprocessingState* preprocessing_state = transformation->preprocessing_state;

  int num_fixed_vars;
  int* fixed_var_indices;
  double* fixed_var_values;

  SLEQP_CALL(sleqp_preprocessing_state_fixed_variables(preprocessing_state,
                                                       &num_fixed_vars,
                                                       &fixed_var_indices,
                                                       &fixed_var_values));

  SLEQP_CALL(sleqp_sparse_vector_to_raw(sleqp_problem_var_ub(problem),
                                        transformation->dense_cache));

  SleqpConvertedBound* converted_bounds;
  int num_converted_bounds;

  SLEQP_CALL(sleqp_preprocessing_state_converted_bounds(preprocessing_state,
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
                                                num_fixed_vars,
                                                fixed_var_indices));

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE transform_linear_constraints(SleqpTransformation* transformation)
{
  SleqpProblem* problem = transformation->original_problem;

  int num_fixed_vars;
  int* fixed_var_indices;
  double* fixed_var_values;

  const int num_linear = sleqp_problem_num_linear_constraints(problem);

  const double zero_eps = sleqp_params_get(transformation->params,
                                           SLEQP_PARAM_ZERO_EPS);

  SleqpPreprocessingState* preprocessing_state = transformation->preprocessing_state;

  SLEQP_CALL(sleqp_preprocessing_state_fixed_variables(preprocessing_state,
                                                       &num_fixed_vars,
                                                       &fixed_var_indices,
                                                       &fixed_var_values));

  SleqpSparseVec* fixed_linear_lb = sleqp_problem_linear_lb(problem);
  SleqpSparseVec* fixed_linear_ub = sleqp_problem_linear_ub(problem);

  if(num_fixed_vars > 0)
  {
    SLEQP_CALL(sleqp_sparse_vector_clear(transformation->fixed_variable_values));

    SLEQP_CALL(sleqp_sparse_vector_reserve(transformation->fixed_variable_values,
                                           num_fixed_vars));

    for(int j = 0; j < num_fixed_vars; ++j)
    {
      SLEQP_CALL(sleqp_sparse_vector_push(transformation->fixed_variable_values,
                                          fixed_var_indices[j],
                                          fixed_var_values[j]));
    }

    SLEQP_CALL(sleqp_sparse_matrix_vector_product(sleqp_problem_linear_coeffs(problem),
                                                  transformation->fixed_variable_values,
                                                  transformation->dense_cache));

    SLEQP_CALL(sleqp_sparse_vector_from_raw(transformation->sparse_cache,
                                            transformation->dense_cache,
                                            num_linear,
                                            zero_eps));

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(sleqp_problem_linear_lb(problem),
                                              transformation->sparse_cache,
                                              1.,
                                              -1.,
                                              zero_eps,
                                              transformation->fixed_linear_lb));

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(sleqp_problem_linear_ub(problem),
                                              transformation->sparse_cache,
                                              1.,
                                              -1.,
                                              zero_eps,
                                              transformation->fixed_linear_ub));

    fixed_linear_lb = transformation->fixed_linear_lb;
    fixed_linear_ub = transformation->fixed_linear_ub;
  }

  int num_removed_cons;
  int* removed_cons_indices;

  SLEQP_CALL(sleqp_preprocessing_state_removed_linear_constraints(preprocessing_state,
                                                                  &num_removed_cons,
                                                                  &removed_cons_indices));

  SLEQP_CALL(sleqp_preprocessing_remove_entries(fixed_linear_lb,
                                                transformation->transformed_linear_lb,
                                                num_removed_cons,
                                                removed_cons_indices));

  SLEQP_CALL(sleqp_preprocessing_remove_entries(fixed_linear_ub,
                                                transformation->transformed_linear_ub,
                                                num_removed_cons,
                                                removed_cons_indices));

  SLEQP_CALL(sleqp_preprocessing_remove_matrix_entries(sleqp_problem_linear_coeffs(problem),
                                                       transformation->transformed_linear_coeffs,
                                                       num_fixed_vars,
                                                       fixed_var_indices,
                                                       num_removed_cons,
                                                       removed_cons_indices));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_transformation_create_transformed_problem(SleqpTransformation* transformation,
                                                              SleqpProblem** star)
{
  SleqpProblem* problem = transformation->original_problem;

  SleqpPreprocessingState* preprocessing_state = transformation->preprocessing_state;

  const int num_fixed_vars = sleqp_preprocessing_state_num_fixed_variables(preprocessing_state);

  const int num_removed_cons = sleqp_preprocessing_state_num_removed_linear_constraints(preprocessing_state);

  SLEQP_CALL(create_transformed_func(transformation,
                                     &transformation->transformed_func));

  SLEQP_CALL(create_transformed_var_lb(transformation));

  SLEQP_CALL(create_transformed_var_ub(transformation));

  SLEQP_CALL(transform_linear_constraints(transformation));

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

    sleqp_assert(num_transformed_variables + num_fixed_vars == num_variables);
    sleqp_assert(num_transformed_constraints + num_removed_cons == num_constraints);
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

  SLEQP_CALL(sleqp_sparse_vector_free(&transformation->fixed_linear_ub));
  SLEQP_CALL(sleqp_sparse_vector_free(&transformation->fixed_linear_lb));

  SLEQP_CALL(sleqp_sparse_vector_free(&transformation->fixed_variable_values));

  SLEQP_CALL(sleqp_sparse_vector_free(&transformation->transformed_var_ub));
  SLEQP_CALL(sleqp_sparse_vector_free(&transformation->transformed_var_lb));

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
