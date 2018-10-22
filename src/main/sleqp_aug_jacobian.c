#include "sleqp_aug_jacobian.h"

#include "sleqp_mem.h"

#include "sparse/sleqp_sparse_factorization.h"

struct SleqpAugJacobian
{
  SleqpProblem* problem;

  SleqpSparseMatrix* augmented_matrix;
  SleqpSparseFactorization* factorization;

  int num_active_conss;
  int num_active_vars;

  int active_set_size;

  // a mapping of 0..num_variables - 1 -> pos in the active set or -1
  int* active_vars;

  // a mapping of 0..num_constraints - 1 -> pos in the active set or -1
  int* active_conss;

  int* set_indices;

  int* col_indices;

};

static SLEQP_RETCODE fill_augmented_jacobian(SleqpAugJacobian* jacobian,
                                             SleqpIterate* iterate,
                                             int total_nnz)
{
  SleqpProblem* problem = jacobian->problem;

  int num_active_vars = jacobian->num_active_vars;
  int active_set_size = jacobian->active_set_size;

  SleqpSparseMatrix* cons_jac = iterate->cons_jac;

  SleqpSparseMatrix* augmented_matrix = jacobian->augmented_matrix;

  int num_variables = problem->num_variables;

  SLEQP_ACTIVE_STATE* cons_states = sleqp_active_set_cons_states(iterate->active_set);
  SLEQP_ACTIVE_STATE* var_states = sleqp_active_set_var_states(iterate->active_set);

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
      // push the jacobian of the variable, if it active (given by +/- e_j)

      int variable_index = sleqp_aug_jacobian_variable_index(jacobian, column);

      if(variable_index != -1)
      {
        assert(var_states[column] != SLEQP_INACTIVE);

        bool is_upper = (var_states[column] == SLEQP_ACTIVE_UPPER);

        SLEQP_CALL(sleqp_sparse_matrix_push(augmented_matrix,
                                            num_variables + variable_index,
                                            column,
                                            is_upper ? 1. : -1.));
      }
      else
      {
        assert(var_states[column] == SLEQP_INACTIVE);
      }
    }

    for(int index = cons_jac->cols[column];
        index < cons_jac->cols[column + 1];
        ++index)
    {
      int jac_row = cons_jac->rows[index];
      int act_cons_index = sleqp_aug_jacobian_constraint_index(jacobian, jac_row);

      if(act_cons_index != -1)
      {
        assert(cons_states[jac_row] != SLEQP_INACTIVE);

        bool is_upper = (cons_states[jac_row] == SLEQP_ACTIVE_UPPER);

        double value = is_upper ? cons_jac->data[index] : -(cons_jac->data[index]);

        int aug_row = num_variables + num_active_vars + act_cons_index;

        SLEQP_CALL(sleqp_sparse_matrix_push(augmented_matrix,
                                            aug_row, // row
                                            column, // column
                                            value));

      }
      else
      {
        assert(cons_states[jac_row] == SLEQP_INACTIVE);
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
                                        SleqpProblem* problem)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpAugJacobian* jacobian = *star;

  jacobian->problem = problem;

  SLEQP_CALL(sleqp_sparse_matrix_create(&jacobian->augmented_matrix, 0, 0, 0));

  jacobian->factorization = NULL;

  jacobian->active_set_size = 0;

  SLEQP_CALL(sleqp_calloc(&jacobian->active_vars, problem->num_variables));
  SLEQP_CALL(sleqp_calloc(&jacobian->active_conss, problem->num_constraints));

  int max_set_size = problem->num_constraints + problem->num_variables;
  int max_num_cols = problem->num_variables + max_set_size;

  SLEQP_CALL(sleqp_calloc(&jacobian->set_indices, max_set_size));
  SLEQP_CALL(sleqp_calloc(&jacobian->col_indices, max_num_cols + 1));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_aug_jacobian_free(SleqpAugJacobian** star)
{
  SleqpAugJacobian* jacobian = *star;

  sleqp_free(&jacobian->col_indices);
  sleqp_free(&jacobian->set_indices);

  sleqp_free(&jacobian->active_conss);
  sleqp_free(&jacobian->active_vars);

  if(jacobian->factorization)
  {
    SLEQP_CALL(sleqp_sparse_factorization_free(&jacobian->factorization));
  }

  sleqp_sparse_matrix_free(&jacobian->augmented_matrix);

  sleqp_free(star);

  return SLEQP_OKAY;
}

int sleqp_aug_jacobian_num_active_vars(SleqpAugJacobian* jacobian)
{
  return jacobian->num_active_vars;
}

int sleqp_aug_jacobian_num_active_conss(SleqpAugJacobian* jacobian)
{
  return jacobian->num_active_conss;
}

int sleqp_aug_jacobian_active_set_size(SleqpAugJacobian* jacobian)
{
  return jacobian->active_set_size;
}

SLEQP_RETCODE sleqp_aug_jacobian_set_iterate(SleqpAugJacobian* jacobian,
                                             SleqpIterate* iterate)
{
  SleqpProblem* problem = jacobian->problem;

  SLEQP_ACTIVE_STATE* cons_states = sleqp_active_set_cons_states(iterate->active_set);
  SLEQP_ACTIVE_STATE* var_states = sleqp_active_set_var_states(iterate->active_set);

  jacobian->num_active_conss = 0;
  jacobian->num_active_vars = 0;
  jacobian->active_set_size = 0;

  // first step: count active vars / cons
  for(int j = 0; j < problem->num_variables; ++j)
  {
    if(var_states[j] != SLEQP_INACTIVE)
    {
      jacobian->active_vars[j] = (jacobian->num_active_vars)++;

      jacobian->set_indices[(jacobian->active_set_size)++] = j;
    }
    else
    {
      jacobian->active_vars[j] = -1;
    }
  }

  for(int i = 0; i < problem->num_constraints; ++i)
  {
    if(cons_states[i] != SLEQP_INACTIVE)
    {
      jacobian->active_conss[i] = jacobian->num_active_conss++;

      jacobian->set_indices[jacobian->active_set_size++] = problem->num_variables + i;
    }
    else
    {
      jacobian->active_conss[i] = -1;
    }
  }

  jacobian->active_set_size = jacobian->num_active_conss + jacobian->num_active_vars;

  // we overestimate here...
  int constraint_nnz = iterate->cons_jac->nnz;

  int variable_nnz = jacobian->num_active_vars;

  int total_nnz = problem->num_variables // identity
    + 2*(constraint_nnz + variable_nnz); // cons jac nnz + active variables

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

int sleqp_aug_jacobian_constraint_index(SleqpAugJacobian* jacobian,
                                        int index)
{
  return jacobian->active_conss[index];
}

int sleqp_aug_jacobian_variable_index(SleqpAugJacobian* jacobian,
                                      int index)
{
  return jacobian->active_vars[index];
}

int sleqp_aug_jacobian_get_set_index(SleqpAugJacobian* jacobian,
                                     int index)
{
  return jacobian->set_indices[index];
}

int sleqp_aug_jacobian_size(SleqpAugJacobian* jacobian)
{
  return jacobian->problem->num_variables + jacobian->active_set_size;
}

SLEQP_RETCODE sleqp_aug_jacobian_min_norm_solution(SleqpAugJacobian* jacobian,
                                                   SleqpSparseVec* rhs,
                                                   SleqpSparseVec* sol)
{
  assert(jacobian->factorization);

  SleqpProblem* problem = jacobian->problem;
  SleqpSparseFactorization* factorization = jacobian->factorization;

  int num_variables = problem->num_variables;

  assert(sol->dim == num_variables);
  assert(rhs->dim == jacobian->active_set_size);

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
                                                problem->num_variables));

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
                                                  num_variables));
  }

  if(dual_sol)
  {
    assert(dual_sol->dim == active_set_size);

    SLEQP_CALL(sleqp_sparse_factorization_get_sol(factorization,
                                                  dual_sol,
                                                  num_variables,
                                                  total_size));

    //SLEQP_CALL(sleqp_sparse_vector_scale(dual_sol, -1.));

  }

  // erase the zeros
  SLEQP_CALL(sleqp_sparse_vector_resize(rhs, num_variables));

  return SLEQP_OKAY;
}
