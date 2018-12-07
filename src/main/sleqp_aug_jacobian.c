#include "sleqp_aug_jacobian.h"

#include "sleqp_mem.h"

#include "sparse/sleqp_sparse_factorization.h"

struct SleqpAugJacobian
{
  SleqpProblem* problem;
  SleqpParams* params;

  int active_set_size;
  int max_set_size;

  SleqpSparseMatrix* augmented_matrix;
  SleqpSparseFactorization* factorization;

  int* col_indices;
};

static SLEQP_RETCODE fill_augmented_jacobian(SleqpAugJacobian* jacobian,
                                             SleqpIterate* iterate,
                                             int total_nnz)
{
  SleqpProblem* problem = jacobian->problem;

  SleqpSparseMatrix* cons_jac = iterate->cons_jac;

  SleqpSparseMatrix* augmented_matrix = jacobian->augmented_matrix;

  SleqpActiveSet* active_set = iterate->active_set;

  const int num_variables = problem->num_variables;
  const int active_set_size = sleqp_active_set_size(active_set);

  int augmented_size = num_variables + active_set_size;

  SLEQP_CALL(sleqp_sparse_matrix_resize(augmented_matrix,
                                        augmented_size,
                                        augmented_size));

  for(int column = 0; column < num_variables; ++column)
  {
    SLEQP_CALL(sleqp_sparse_matrix_add_column(augmented_matrix,
                                              column));

    // push identity part first...
    SLEQP_CALL(sleqp_sparse_matrix_push(augmented_matrix,
                                        column,
                                        column,
                                        1.));

    {
      // push the jacobian of the variable, if active

      int variable_index = sleqp_active_set_get_variable_index(active_set, column);

      SLEQP_ACTIVE_STATE var_state = sleqp_active_set_get_variable_state(active_set, column);

      if(variable_index != -1)
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

    for(int index = cons_jac->cols[column];
        index < cons_jac->cols[column + 1];
        ++index)
    {
      int jac_row = cons_jac->rows[index];

      int act_cons_index = sleqp_active_set_get_constraint_index(active_set, jac_row);

      SLEQP_ACTIVE_STATE cons_state = sleqp_active_set_get_constraint_state(active_set, jac_row);

      if(act_cons_index != -1)
      {
        assert(cons_state != SLEQP_INACTIVE);

        const double value = cons_jac->data[index];

        int aug_row = num_variables + act_cons_index;

        SLEQP_CALL(sleqp_sparse_matrix_push(augmented_matrix,
                                            aug_row, // row
                                            column, // column
                                            value));

      }
      else
      {
        assert(cons_state == SLEQP_INACTIVE);
      }
    }
  }

  for(int column = num_variables + 1; column < augmented_matrix->num_cols + 1; ++column)
  {
    augmented_matrix->cols[column] = 0;
  }

  // count the elements in the columns
  for(int column = 0; column < num_variables; ++column)
  {
    for(int index = augmented_matrix->cols[column]; index < augmented_matrix->cols[column + 1]; ++index)
    {
      if(augmented_matrix->rows[index] < num_variables)
      {
        continue;
      }

      ++augmented_matrix->cols[augmented_matrix->rows[index] + 1];
    }
  }

  // construct the cumulative sum
  for(int column = num_variables + 1; column < augmented_matrix->num_cols + 1; ++column)
  {
    augmented_matrix->cols[column] += augmented_matrix->cols[column - 1];
    jacobian->col_indices[column] = 0;
  }

  // insert the entries
  for(int column = 0; column < num_variables; ++column)
  {
    for(int index = augmented_matrix->cols[column]; index < augmented_matrix->cols[column + 1]; ++index)
    {
      if(augmented_matrix->rows[index] < num_variables)
      {
        continue;
      }

      // target column : augmented_matrix->rows[index]
      // target row: *column*
      // target data: augmented_matrix->data[index]
      // target index: augmented_matrix->cols[*target column*] + (col_indices[*target column*]++)

      int target_column = augmented_matrix->rows[index];
      int target_index = augmented_matrix->cols[target_column] + jacobian->col_indices[target_column + 1]++;

      augmented_matrix->data[target_index] = augmented_matrix->data[index];
      augmented_matrix->rows[target_index] = column;
      ++augmented_matrix->nnz;
    }
  }

  assert(sleqp_sparse_matrix_valid(augmented_matrix));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_aug_jacobian_create(SleqpAugJacobian** star,
                                        SleqpProblem* problem,
                                        SleqpParams* params)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpAugJacobian* jacobian = *star;

  jacobian->problem = problem;
  jacobian->params = params;

  jacobian->active_set_size = 0;

  SLEQP_CALL(sleqp_sparse_matrix_create(&jacobian->augmented_matrix, 0, 0, 0));

  jacobian->factorization = NULL;

  jacobian->max_set_size = problem->num_constraints + problem->num_variables;
  int max_num_cols = problem->num_variables + jacobian->max_set_size;

  SLEQP_CALL(sleqp_calloc(&jacobian->col_indices, max_num_cols + 1));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_aug_jacobian_free(SleqpAugJacobian** star)
{
  SleqpAugJacobian* jacobian = *star;

  sleqp_free(&jacobian->col_indices);

  if(jacobian->factorization)
  {
    SLEQP_CALL(sleqp_sparse_factorization_free(&jacobian->factorization));
  }

  sleqp_sparse_matrix_free(&jacobian->augmented_matrix);

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_aug_jacobian_set_iterate(SleqpAugJacobian* jacobian,
                                             SleqpIterate* iterate)
{
  SleqpProblem* problem = jacobian->problem;
  SleqpActiveSet* active_set = iterate->active_set;

  jacobian->active_set_size = sleqp_active_set_size(active_set);

  // we overestimate here...
  int constraint_nnz = iterate->cons_jac->nnz;

  int variable_nnz = sleqp_active_set_num_active_vars(active_set);

  int total_nnz = problem->num_variables // identity
    + 2*(constraint_nnz + variable_nnz); // cons jac nnz + active variables

  SLEQP_CALL(sleqp_sparse_matrix_clear(jacobian->augmented_matrix));

  SLEQP_CALL(sleqp_sparse_matrix_reserve(jacobian->augmented_matrix, total_nnz));

  SLEQP_CALL(fill_augmented_jacobian(jacobian,
                                     iterate,
                                     total_nnz));

  if(jacobian->factorization)
  {
    SLEQP_CALL(sleqp_sparse_factorization_free(&jacobian->factorization));
  }

  SLEQP_CALL(sleqp_sparse_factorization_create(&jacobian->factorization,
                                               jacobian->augmented_matrix));

  {
    SleqpSparseMatrix* matrix = jacobian->augmented_matrix;

    int num_variables = problem->num_variables;
    int total_size = num_variables + jacobian->active_set_size;

    assert(matrix->num_cols == matrix->num_rows);
    assert(matrix->num_cols == total_size);
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_aug_jacobian_min_norm_solution(SleqpAugJacobian* jacobian,
                                                   SleqpSparseVec* rhs,
                                                   SleqpSparseVec* sol)
{
  assert(jacobian->factorization);

  SleqpProblem* problem = jacobian->problem;
  SleqpSparseFactorization* factorization = jacobian->factorization;

  double zero_eps = sleqp_params_get_zero_eps(jacobian->params);

  int num_variables = problem->num_variables;

  assert(sol->dim == num_variables);

  rhs->dim += num_variables;

  for(int k = 0; k < rhs->nnz; ++k)
  {
    rhs->indices[k] += num_variables;
  }

  SLEQP_CALL(sleqp_sparse_factorization_solve(factorization,
                                              rhs));

  SLEQP_CALL(sleqp_sparse_factorization_get_sol(factorization,
                                                sol,
                                                0,
                                                problem->num_variables,
                                                zero_eps));

  rhs->dim -= num_variables;

  for(int k = 0; k < rhs->nnz; ++k)
  {
    rhs->indices[k] -= num_variables;
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_aug_jacobian_projection(SleqpAugJacobian* jacobian,
                                            SleqpSparseVec* rhs,
                                            SleqpSparseVec* primal_sol,
                                            SleqpSparseVec* dual_sol)
{
  assert(jacobian->factorization);

  SleqpProblem* problem = jacobian->problem;
  SleqpSparseFactorization* factorization = jacobian->factorization;

  double zero_eps = sleqp_params_get_zero_eps(jacobian->params);

  int num_variables = problem->num_variables;
  int active_set_size = jacobian->active_set_size;
  int total_size = num_variables + active_set_size;

  assert(rhs->dim == num_variables);

  // just add some zeros...
  SLEQP_CALL(sleqp_sparse_vector_resize(rhs, total_size));

  SLEQP_CALL(sleqp_sparse_factorization_solve(factorization,
                                              rhs));

  if(primal_sol)
  {
    assert(primal_sol->dim == num_variables);

    SLEQP_CALL(sleqp_sparse_factorization_get_sol(factorization,
                                                  primal_sol,
                                                  0,
                                                  num_variables,
                                                  zero_eps));
  }

  if(dual_sol)
  {
    assert(dual_sol->dim == active_set_size);

    SLEQP_CALL(sleqp_sparse_factorization_get_sol(factorization,
                                                  dual_sol,
                                                  num_variables,
                                                  total_size,
                                                  zero_eps));

    //SLEQP_CALL(sleqp_sparse_vector_scale(dual_sol, -1.));

  }

  // erase the zeros
  SLEQP_CALL(sleqp_sparse_vector_resize(rhs, num_variables));

  return SLEQP_OKAY;
}
