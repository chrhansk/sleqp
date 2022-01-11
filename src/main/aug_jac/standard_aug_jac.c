#include "standard_aug_jac.h"

#include <assert.h>

#include "log.h"
#include "mem.h"
#include "problem.h"
#include "working_set.h"

#include "factorization/factorization.h"

typedef struct
{
  SleqpProblem* problem;
  SleqpParams* params;

  int working_set_size;
  int max_set_size;

  bool has_factorization;

  SleqpSparseMatrix* augmented_matrix;
  SleqpFactorization* factorization;

  SleqpWorkingSet* working_set;

  SleqpTimer* factorization_timer;
  SleqpTimer* substitution_timer;

  double condition_estimate;

  int* col_indices;
} AugJacData;

static SLEQP_RETCODE
fill_augmented_jacobian(AugJacData* jacobian,
                        SleqpIterate* iterate,
                        int total_nnz)
{
  SleqpProblem* problem = jacobian->problem;

  SleqpSparseMatrix* cons_jac = sleqp_iterate_cons_jac(iterate);

  SleqpSparseMatrix* augmented_matrix = jacobian->augmented_matrix;

  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);

  const int num_variables    = sleqp_problem_num_vars(problem);
  const int working_set_size = sleqp_working_set_size(working_set);

  int augmented_size = num_variables + working_set_size;

  SLEQP_CALL(sleqp_sparse_matrix_resize(augmented_matrix,
                                        augmented_size,
                                        augmented_size));

  for (int column = 0; column < num_variables; ++column)
  {
    SLEQP_CALL(sleqp_sparse_matrix_push_column(augmented_matrix, column));

    // push identity part first...
    SLEQP_CALL(sleqp_sparse_matrix_push(augmented_matrix, column, column, 1.));

    {
      // push the jacobian of the variable, if active

      int variable_index = sleqp_working_set_var_index(working_set, column);

      SLEQP_ACTIVE_STATE var_state
        = sleqp_working_set_var_state(working_set, column);

      if (variable_index != -1)
      {
        assert(var_state != SLEQP_INACTIVE);

        SLEQP_CALL(sleqp_sparse_matrix_push(augmented_matrix,
                                            num_variables + variable_index,
                                            column,
                                            1.));
      }
      else
      {
        assert(var_state == SLEQP_INACTIVE);
      }
    }

    int* cons_jac_cols    = sleqp_sparse_matrix_cols(cons_jac);
    int* cons_jac_rows    = sleqp_sparse_matrix_rows(cons_jac);
    double* cons_jac_data = sleqp_sparse_matrix_data(cons_jac);

    for (int index = cons_jac_cols[column]; index < cons_jac_cols[column + 1];
         ++index)
    {
      int jac_row = cons_jac_rows[index];

      int act_cons_index = sleqp_working_set_cons_index(working_set, jac_row);

      SLEQP_ACTIVE_STATE cons_state
        = sleqp_working_set_cons_state(working_set, jac_row);

      if (act_cons_index != -1)
      {
        assert(cons_state != SLEQP_INACTIVE);

        const double value = cons_jac_data[index];

        int aug_row = num_variables + act_cons_index;

        SLEQP_CALL(sleqp_sparse_matrix_push(augmented_matrix,
                                            aug_row, // row
                                            column,  // column
                                            value));
      }
      else
      {
        assert(cons_state == SLEQP_INACTIVE);
      }
    }
  }

  int aug_num_cols = sleqp_sparse_matrix_num_cols(augmented_matrix);
  int* aug_cols    = sleqp_sparse_matrix_cols(augmented_matrix);
  int* aug_rows    = sleqp_sparse_matrix_rows(augmented_matrix);
  double* aug_data = sleqp_sparse_matrix_data(augmented_matrix);

  for (int column = num_variables + 1; column < aug_num_cols + 1; ++column)
  {
    aug_cols[column] = 0;
  }

  // count the elements in the columns
  for (int column = 0; column < num_variables; ++column)
  {
    for (int index = aug_cols[column]; index < aug_cols[column + 1]; ++index)
    {
      if (aug_rows[index] < num_variables)
      {
        continue;
      }

      ++aug_cols[aug_rows[index] + 1];
    }
  }

  // construct the cumulative sum
  for (int column = num_variables + 1; column < aug_num_cols + 1; ++column)
  {
    aug_cols[column] += aug_cols[column - 1];
    jacobian->col_indices[column] = 0;
  }

  int aug_total_nnz = sleqp_sparse_matrix_nnz(augmented_matrix);

  // insert the entries
  for (int column = 0; column < num_variables; ++column)
  {
    for (int index = aug_cols[column]; index < aug_cols[column + 1]; ++index)
    {
      if (aug_rows[index] < num_variables)
      {
        continue;
      }

      // target column : aug_rows[index]
      // target row: *column*
      // target data: augmented_matrix->data[index]
      // target index: aug_cols[*target column*] + (col_indices[*target
      // column*]++)

      int target_column = aug_rows[index];
      int target_index
        = aug_cols[target_column] + jacobian->col_indices[target_column + 1]++;

      aug_data[target_index] = aug_data[index];
      aug_rows[target_index] = column;
      ++aug_total_nnz;
    }
  }

  SLEQP_CALL(sleqp_sparse_matrix_set_nnz(augmented_matrix, aug_total_nnz));

  assert(sleqp_sparse_matrix_is_valid(augmented_matrix));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
aug_jac_set_iterate(SleqpIterate* iterate, void* data)
{
  AugJacData* jacobian = (AugJacData*)data;

  SleqpProblem* problem        = jacobian->problem;
  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);

  if (jacobian->working_set)
  {
    if (sleqp_working_set_eq(working_set, jacobian->working_set)
        && jacobian->has_factorization)
    {
      return SLEQP_OKAY;
    }
    else
    {
      SLEQP_CALL(sleqp_working_set_copy(working_set, jacobian->working_set));
    }
  }

  const int num_variables = sleqp_problem_num_vars(problem);

  jacobian->working_set_size = sleqp_working_set_size(working_set);

  jacobian->condition_estimate = SLEQP_NONE;

  // we overestimate here...
  int constraint_nnz = sleqp_sparse_matrix_nnz(sleqp_iterate_cons_jac(iterate));

  int variable_nnz = sleqp_working_set_num_active_vars(working_set);

  int total_nnz
    = num_variables                          // identity
      + 2 * (constraint_nnz + variable_nnz); // cons jac nnz + active variables

  SLEQP_CALL(sleqp_sparse_matrix_clear(jacobian->augmented_matrix));

  SLEQP_CALL(
    sleqp_sparse_matrix_reserve(jacobian->augmented_matrix, total_nnz));

  SLEQP_CALL(fill_augmented_jacobian(jacobian, iterate, total_nnz));

  SLEQP_CALL(sleqp_timer_start(jacobian->factorization_timer));

  SLEQP_CALL(sleqp_factorization_set_matrix(jacobian->factorization,
                                            jacobian->augmented_matrix));

  jacobian->has_factorization = true;

  SLEQP_CALL(sleqp_timer_stop(jacobian->factorization_timer));

  SLEQP_CALL(
    sleqp_factorization_condition_estimate(jacobian->factorization,
                                           &jacobian->condition_estimate));

  {
    SleqpSparseMatrix* matrix = jacobian->augmented_matrix;

    int total_size = num_variables + jacobian->working_set_size;

    assert(sleqp_sparse_matrix_is_quadratic(matrix));
    assert(sleqp_sparse_matrix_num_cols(matrix) == total_size);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
aug_jac_condition(bool* exact, double* condition, void* data)
{
  AugJacData* jacobian = (AugJacData*)data;

  *exact     = false;
  *condition = jacobian->condition_estimate;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
aug_jac_min_norm_solution(const SleqpSparseVec* _rhs,
                          SleqpSparseVec* sol,
                          void* data)
{
  AugJacData* jacobian = (AugJacData*)data;

  // Cast away constness
  SleqpSparseVec* rhs = (SleqpSparseVec*)_rhs;

  assert(jacobian->factorization);

  SLEQP_CALL(sleqp_timer_start(jacobian->substitution_timer));

  SleqpProblem* problem             = jacobian->problem;
  SleqpFactorization* factorization = jacobian->factorization;

  double zero_eps = sleqp_params_value(jacobian->params, SLEQP_PARAM_ZERO_EPS);

  const int num_variables = sleqp_problem_num_vars(problem);

  assert(sol->dim == num_variables);

  rhs->dim += num_variables;

  for (int k = 0; k < rhs->nnz; ++k)
  {
    rhs->indices[k] += num_variables;
  }

  SLEQP_CALL(sleqp_factorization_solve(factorization, rhs));

  SLEQP_CALL(sleqp_factorization_solution(factorization,
                                          sol,
                                          0,
                                          num_variables,
                                          zero_eps));

  rhs->dim -= num_variables;

  for (int k = 0; k < rhs->nnz; ++k)
  {
    rhs->indices[k] -= num_variables;
  }

  SLEQP_CALL(sleqp_timer_stop(jacobian->substitution_timer));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
aug_jac_projection(const SleqpSparseVec* _rhs,
                   SleqpSparseVec* primal_sol,
                   SleqpSparseVec* dual_sol,
                   void* data)
{
  AugJacData* jacobian = (AugJacData*)data;

  // Cast away constness
  SleqpSparseVec* rhs = (SleqpSparseVec*)_rhs;

  assert(jacobian->factorization);

  SLEQP_CALL(sleqp_timer_start(jacobian->substitution_timer));

  SleqpProblem* problem             = jacobian->problem;
  SleqpFactorization* factorization = jacobian->factorization;

  double zero_eps = sleqp_params_value(jacobian->params, SLEQP_PARAM_ZERO_EPS);

  const int num_variables    = sleqp_problem_num_vars(problem);
  const int working_set_size = jacobian->working_set_size;
  const int total_size       = num_variables + working_set_size;

  assert(rhs->dim == num_variables);

  // just add some zeros...
  SLEQP_CALL(sleqp_sparse_vector_resize(rhs, total_size));

  SLEQP_CALL(sleqp_factorization_solve(factorization, rhs));

  if (primal_sol)
  {
    assert(primal_sol->dim == num_variables);

    SLEQP_CALL(sleqp_factorization_solution(factorization,
                                            primal_sol,
                                            0,
                                            num_variables,
                                            zero_eps));
  }

  if (dual_sol)
  {
    assert(dual_sol->dim == working_set_size);

    SLEQP_CALL(sleqp_factorization_solution(factorization,
                                            dual_sol,
                                            num_variables,
                                            total_size,
                                            zero_eps));

    // SLEQP_CALL(sleqp_sparse_vector_scale(dual_sol, -1.));
  }

  // erase the zeros
  SLEQP_CALL(sleqp_sparse_vector_resize(rhs, num_variables));

  SLEQP_CALL(sleqp_timer_stop(jacobian->substitution_timer));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
aug_jac_free(void* data)
{
  AugJacData* jacobian = (AugJacData*)data;

  SLEQP_CALL(sleqp_timer_free(&jacobian->substitution_timer));

  SLEQP_CALL(sleqp_timer_free(&jacobian->factorization_timer));

  sleqp_free(&jacobian->col_indices);

  SLEQP_CALL(sleqp_working_set_release(&jacobian->working_set));

  SLEQP_CALL(sleqp_factorization_release(&jacobian->factorization));

  SLEQP_CALL(sleqp_sparse_matrix_release(&jacobian->augmented_matrix));

  SLEQP_CALL(sleqp_params_release(&jacobian->params));

  SLEQP_CALL(sleqp_problem_release(&jacobian->problem));

  sleqp_free(&jacobian);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
aug_jacobian_data_create(AugJacData** star,
                         SleqpProblem* problem,
                         SleqpParams* params,
                         SleqpFactorization* factorization)
{
  SLEQP_CALL(sleqp_malloc(star));

  const int num_variables   = sleqp_problem_num_vars(problem);
  const int num_constraints = sleqp_problem_num_cons(problem);

  AugJacData* jacobian = *star;

  *jacobian = (AugJacData){0};

  jacobian->problem = problem;

  SLEQP_CALL(sleqp_problem_capture(jacobian->problem));

  SLEQP_CALL(sleqp_params_capture(params));

  jacobian->params = params;

  jacobian->condition_estimate = -1;

  SLEQP_CALL(sleqp_sparse_matrix_create(&jacobian->augmented_matrix, 0, 0, 0));

  jacobian->max_set_size = num_constraints + num_variables;
  int max_num_cols       = num_variables + jacobian->max_set_size;

  jacobian->has_factorization = false;

  SLEQP_CALL(sleqp_factorization_capture(factorization));

  jacobian->factorization = factorization;

  const bool fixed_jacobian = !(sleqp_problem_has_nonlinear_cons(problem));

  if (fixed_jacobian)
  {
    SLEQP_CALL(sleqp_working_set_create(&jacobian->working_set, problem));
  }

  SLEQP_CALL(sleqp_alloc_array(&jacobian->col_indices, max_num_cols + 1));

  SLEQP_CALL(sleqp_timer_create(&jacobian->factorization_timer));

  SLEQP_CALL(sleqp_timer_create(&jacobian->substitution_timer));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_standard_aug_jac_create(SleqpAugJac** star,
                              SleqpProblem* problem,
                              SleqpParams* params,
                              SleqpFactorization* factorization)
{
  AugJacData* aug_jac_data;

  SLEQP_CALL(
    aug_jacobian_data_create(&aug_jac_data, problem, params, factorization));

  SleqpAugJacCallbacks callbacks
    = {.set_iterate       = aug_jac_set_iterate,
       .min_norm_solution = aug_jac_min_norm_solution,
       .projection        = aug_jac_projection,
       .condition         = aug_jac_condition,
       .free              = aug_jac_free};

  SLEQP_CALL(sleqp_aug_jac_create(star, problem, &callbacks, aug_jac_data));

  return SLEQP_OKAY;
}
