#include "sleqp_dual_estimation.h"

#include <umfpack.h>

#include "sleqp_cmp.h"
#include "sleqp_mem.h"

struct SleqpDualEstimationData
{
  int num_variables;
  int num_constraints;

  int* active_vars;
  int* active_cons;

  int* col_indices;
  SleqpSparseMatrix* estimation_matrix;
  double* estimation_rhs;
  double* estimation_values;
};


SLEQP_RETCODE sleqp_dual_estimation_data_create(SleqpDualEstimationData** star,
                                                SleqpProblem* problem)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpDualEstimationData* estimation_data = *star;

  estimation_data->num_variables = problem->num_variables;
  estimation_data->num_constraints = problem->num_constraints;

  SLEQP_CALL(sleqp_calloc(&estimation_data->active_vars, estimation_data->num_variables));
  SLEQP_CALL(sleqp_calloc(&estimation_data->active_cons, estimation_data->num_constraints));

  SLEQP_CALL(sleqp_calloc(&estimation_data->col_indices,
                          estimation_data->num_variables + estimation_data->num_constraints));

  SLEQP_CALL(sleqp_sparse_matrix_create(&estimation_data->estimation_matrix,
                                        0,
                                        0,
                                        0));

  int estimation_size_bound = 2*problem->num_variables + problem->num_constraints;

  SLEQP_CALL(sleqp_calloc(&estimation_data->estimation_rhs, estimation_size_bound));
  SLEQP_CALL(sleqp_calloc(&estimation_data->estimation_values, estimation_size_bound));

  for(int i = 0; i < estimation_size_bound; ++i)
  {
    estimation_data->estimation_rhs[i] = 0;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE resize_estimation_matrix(SleqpDualEstimationData* estimation_data,
                                              int dim,
                                              int total_nnz)
{
  SLEQP_CALL(sleqp_sparse_matrix_reserve(estimation_data->estimation_matrix,
                                         total_nnz));

  SLEQP_CALL(sleqp_sparse_matrix_resize(estimation_data->estimation_matrix,
                                        dim,
                                        dim));

  estimation_data->estimation_matrix->nnz = 0;

  return SLEQP_OKAY;
}

static int count_active_constraint_nnz(SleqpSparseMatrix* cons_jac,
                                          SleqpActiveSet* active_set)
{
  int constraint_nnz = 0;
  SLEQP_ACTIVE_STATE* cons_states = sleqp_active_set_cons_states(active_set);

  int col = 0;

  for(int index = 0; index < cons_jac->nnz; ++index)
  {
    while(col < cons_jac->cols[index])
    {
      ++col;
    }

    if(cons_states[index] != SLEQP_INACTIVE)
    {
      continue;
    }

    ++constraint_nnz;
  }

  return constraint_nnz;
}

static SLEQP_RETCODE fill_estimation_matrix(SleqpDualEstimationData* estimation_data,
                                            SleqpIterate* iterate,
                                            int num_active_vars,
                                            int num_active_cons,
                                            int total_nnz)
{
  SleqpSparseMatrix* cons_jac = iterate->cons_jac;
  SleqpSparseMatrix* matrix = estimation_data->estimation_matrix;

  int num_active = num_active_cons + num_active_vars;
  int num_variables = estimation_data->num_variables;

  SLEQP_ACTIVE_STATE* cons_states = sleqp_active_set_cons_states(iterate->active_set);
  SLEQP_ACTIVE_STATE* var_states = sleqp_active_set_var_states(iterate->active_set);

  for(int column = 0; column < num_variables; ++column)
  {
    SLEQP_CALL(sleqp_sparse_matrix_add_column(matrix,
                                              column));

    // push identity part first...
    SLEQP_CALL(sleqp_sparse_matrix_push(matrix,
                                        column,
                                        column,
                                        1.));

    for(int index = cons_jac->cols[column];
        index < cons_jac->cols[column + 1];
        ++index)
    {
      int cons_row = cons_jac->rows[index];

      if(cons_states[cons_row] == SLEQP_INACTIVE)
      {
        continue;
      }

      int is_upper = (cons_states[cons_row] == SLEQP_ACTIVE_UPPER);

      cons_row = estimation_data->active_cons[cons_row];

      SLEQP_CALL(sleqp_sparse_matrix_push(matrix,
                                          num_variables + cons_row, // row
                                          column, // column
                                          is_upper ? cons_jac->data[index] : -(cons_jac->data[index])));
    }

    if(var_states[column] != SLEQP_INACTIVE)
    {
      int is_upper = (var_states[column] == SLEQP_ACTIVE_UPPER);

      SLEQP_CALL(sleqp_sparse_matrix_push(matrix,
                                          num_variables + num_active_cons + estimation_data->active_vars[column],
                                          column,
                                          is_upper ? 1. : -1.));
    }
  }

  for(int column = num_variables + 1; column < matrix->num_cols + 1; ++column)
  {
    matrix->cols[column] = 0;
  }

  // count the elements in the columns
  for(int column = 0; column < num_variables; ++column)
  {
    for(int index = matrix->cols[column]; index < matrix->cols[column + 1]; ++index)
    {
      if(matrix->rows[index] < num_variables)
      {
        continue;
      }

      ++matrix->cols[matrix->rows[index] + 1];
    }
  }

  // construct the cumulative sum
  for(int column = num_variables + 1; column < matrix->num_cols + 1; ++column)
  {
    matrix->cols[column] += matrix->cols[column - 1];
    estimation_data->col_indices[column] = 0;
  }

  // insert the entries
  for(int column = 0; column < num_variables; ++column)
  {
    for(int index = matrix->cols[column]; index < matrix->cols[column + 1]; ++index)
    {
      if(matrix->rows[index] < num_active)
      {
        continue;
      }

      // target column : matrix->rows[index]
      // target row: *column*
      // target data: matrix->data[index]
      // target index: matrix->cols[*target column*] + (col_indices[*target column*]++)

      int target_column = matrix->rows[index];
      int target_index = matrix->cols[target_column] + estimation_data->col_indices[target_column + 1]++;

      matrix->data[target_index] = matrix->data[index];
      matrix->rows[target_index] = column;
      ++matrix->nnz;
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE construct_estimation_rhs(SleqpDualEstimationData* estimation_data,
                                              SleqpIterate* iterate)
{
  SleqpSparseVec* grad = iterate->func_grad;
  double* rhs = estimation_data->estimation_rhs;

  int index = 0;

  for(int i = 0; i < grad->dim; ++i)
  {
    if(index < grad->nnz && grad->indices[index] == i)
    {
      rhs[i] = -(grad->data[index++]);
    }
    else
    {
      rhs[i] = 0;
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE solve_estimation(SleqpDualEstimationData* estimation_data,
                                      SleqpIterate* iterate,
                                      int num_active_vars,
                                      int num_active_cons)
{
  void* Numeric;
  int status;
  void* Symbolic;

  SleqpSparseMatrix* matrix = estimation_data->estimation_matrix;
  double* solution = estimation_data->estimation_values;
  double* rhs = estimation_data->estimation_rhs;

  assert(matrix->num_cols == matrix->num_rows);

  status = umfpack_di_symbolic(matrix->num_cols,
                               matrix->num_rows,
                               matrix->cols,
                               matrix->rows,
                               matrix->data,
                               &Symbolic,
                               NULL,
                               NULL);

  status = umfpack_di_numeric(matrix->cols,
                              matrix->rows,
                              matrix->data,
                              Symbolic,
                              &Numeric,
                              NULL,
                              NULL);

  umfpack_di_free_symbolic(&Symbolic);

  status = umfpack_di_solve(UMFPACK_A,
                            matrix->cols,
                            matrix->rows,
                            matrix->data,
                            solution,
                            rhs,
                            Numeric,
                            NULL,
                            NULL);

  umfpack_di_free_numeric(&Numeric);

  int num_variables = estimation_data->num_variables;
  int num_constraints = estimation_data->num_constraints;

  SLEQP_ACTIVE_STATE* cons_states = sleqp_active_set_cons_states(iterate->active_set);
  SLEQP_ACTIVE_STATE* var_states = sleqp_active_set_var_states(iterate->active_set);

  SLEQP_CALL(sleqp_sparse_vector_reserve(iterate->cons_dual,
                                         num_active_cons));

  for(int i = 0; i < num_constraints; ++i)
  {
    if(cons_states[i] == SLEQP_ACTIVE_UPPER)
    {
      double value = SLEQP_MAX(solution[num_variables + estimation_data->active_cons[i]], 0);

      if(!sleqp_zero(value))
      {
        SLEQP_CALL(sleqp_sparse_vector_push(iterate->cons_dual,
                                            i,
                                            value));
      }
    }
    else if(cons_states[i] == SLEQP_ACTIVE_LOWER)
    {
      double value = SLEQP_MIN(solution[num_variables + estimation_data->active_cons[i]], 0);

      if(!sleqp_zero(value))
      {
        SLEQP_CALL(sleqp_sparse_vector_push(iterate->cons_dual,
                                            i,
                                            value));
      }
    }
  }

  SLEQP_CALL(sleqp_sparse_vector_reserve(iterate->vars_dual,
                                         num_active_vars));

  for(int i = 0; i < num_variables; ++i)
  {
    if(var_states[i] == SLEQP_ACTIVE_UPPER)
    {
      double value = SLEQP_MAX(-solution[num_variables + num_constraints + estimation_data->active_vars[i]], 0);

      if(!sleqp_zero(value))
      {
        SLEQP_CALL(sleqp_sparse_vector_push(iterate->vars_dual,
                                            i,
                                            value));
      }
    }
    else if(var_states[i] == SLEQP_ACTIVE_LOWER)
    {
      double value = SLEQP_MIN(-solution[num_variables + num_constraints + estimation_data->active_vars[i]], 0);

      if(!sleqp_zero(value))
      {
        SLEQP_CALL(sleqp_sparse_vector_push(iterate->vars_dual,
                                            i,
                                            value));
      }
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_dual_estimation_compute(SleqpDualEstimationData* estimation_data,
                                            SleqpIterate* iterate)
{
  SLEQP_ACTIVE_STATE* cons_states = sleqp_active_set_cons_states(iterate->active_set);
  SLEQP_ACTIVE_STATE* var_states = sleqp_active_set_var_states(iterate->active_set);

  int num_active_cons = 0, num_active_vars = 0;

  // first step: count active cons / vars
  for(int i = 0; i < estimation_data->num_constraints; ++i)
  {
    estimation_data->active_cons[i] = (cons_states[i] != SLEQP_INACTIVE) ? num_active_cons++ : -1;
  }

  for(int i = 0; i < estimation_data->num_variables; ++i)
  {
    estimation_data->active_vars[i] = (var_states[i] != SLEQP_INACTIVE) ? num_active_vars++ : -1;
  }

  int num_active = num_active_vars + num_active_cons;

  int constraint_nnz = count_active_constraint_nnz(iterate->cons_jac,
                                                      iterate->active_set);

  int variable_nnz = num_active_vars;

  int total_nnz = estimation_data->num_variables // identity
    + 2*(constraint_nnz + variable_nnz); // cons jac nnz + active variables

  int matrix_dim = estimation_data->num_variables + num_active;

  SLEQP_CALL(resize_estimation_matrix(estimation_data,
                                      matrix_dim,
                                      total_nnz));

  SLEQP_CALL(fill_estimation_matrix(estimation_data,
                                    iterate,
                                    num_active_vars,
                                    num_active_cons,
                                    total_nnz));

  SLEQP_CALL(construct_estimation_rhs(estimation_data,
                                      iterate));

  SLEQP_CALL(solve_estimation(estimation_data,
                              iterate,
                              num_active_vars,
                              num_active_cons));


  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_dual_estimation_data_free(SleqpDualEstimationData** star)
{
  SleqpDualEstimationData* estimation_data = *star;

  sleqp_free(&estimation_data->estimation_values);
  sleqp_free(&estimation_data->estimation_rhs);

  sleqp_sparse_matrix_free(&estimation_data->estimation_matrix);

  sleqp_free(&estimation_data->col_indices);

  sleqp_free(&estimation_data->active_vars);
  sleqp_free(&estimation_data->active_cons);

  sleqp_free(star);

  return SLEQP_OKAY;
}
