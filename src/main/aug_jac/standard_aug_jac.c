#include "standard_aug_jac.h"

#include "fail.h"
#include "log.h"
#include "mem.h"
#include "problem.h"
#include "working_set.h"

#include "fact/fact.h"

typedef struct
{
  SleqpProblem* problem;
  SleqpParams* params;

  int working_set_size;
  int max_set_size;

  bool has_factorization;

  SleqpMat* augmented_matrix;
  SleqpFact* fact;

  SleqpWorkingSet* working_set;

  SleqpTimer* factorization_timer;
  SleqpTimer* substitution_timer;

  double condition;

  int* col_indices;
} AugJacData;

static SLEQP_RETCODE
add_upper(AugJacData* jacobian)
{
  SleqpProblem* problem      = jacobian->problem;
  SleqpMat* augmented_matrix = jacobian->augmented_matrix;

  const int num_variables = sleqp_problem_num_vars(problem);

  int aug_num_cols = sleqp_mat_num_cols(augmented_matrix);
  int* aug_cols    = sleqp_mat_cols(augmented_matrix);
  int* aug_rows    = sleqp_mat_rows(augmented_matrix);
  double* aug_data = sleqp_mat_data(augmented_matrix);

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

  int aug_total_nnz = sleqp_mat_nnz(augmented_matrix);

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

  SLEQP_CALL(sleqp_mat_set_nnz(augmented_matrix, aug_total_nnz));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
reserve_aug_jac(AugJacData* jacobian, SleqpIterate* iterate, bool lower_only)
{
  SleqpProblem* problem        = jacobian->problem;
  SleqpMat* cons_jac           = sleqp_iterate_cons_jac(iterate);
  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);

  const int num_vars = sleqp_problem_num_vars(problem);
  // we overestimate here...
  const int cons_nnz = sleqp_mat_nnz(cons_jac);

  const int active_var_nnz = sleqp_working_set_num_active_vars(working_set);

  int max_nnz = num_vars + (cons_nnz + active_var_nnz);

  // identity +
  // [2] * (cons jac + active variables)
  if (!lower_only)
  {
    max_nnz += (cons_nnz + active_var_nnz);
  }

  SleqpMat* augmented_matrix = jacobian->augmented_matrix;

  SLEQP_CALL(sleqp_mat_reserve(augmented_matrix, max_nnz));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
fill_aug_jac(AugJacData* jacobian, SleqpIterate* iterate, bool lower_only)
{
  SleqpProblem* problem = jacobian->problem;

  SleqpMat* cons_jac         = sleqp_iterate_cons_jac(iterate);
  SleqpMat* augmented_matrix = jacobian->augmented_matrix;

  SLEQP_CALL(sleqp_mat_clear(augmented_matrix));

  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);

  const int num_variables    = sleqp_problem_num_vars(problem);
  const int working_set_size = sleqp_working_set_size(working_set);

  int augmented_size = num_variables + working_set_size;

  SLEQP_CALL(
    sleqp_mat_resize(augmented_matrix, augmented_size, augmented_size));

  for (int column = 0; column < num_variables; ++column)
  {
    SLEQP_CALL(sleqp_mat_push_col(augmented_matrix, column));

    // push identity part first...
    SLEQP_CALL(sleqp_mat_push(augmented_matrix, column, column, 1.));

    {
      // push the jacobian of the variable, if active

      int variable_index = sleqp_working_set_var_index(working_set, column);

      SLEQP_ACTIVE_STATE var_state
        = sleqp_working_set_var_state(working_set, column);

      if (variable_index != -1)
      {
        assert(var_state != SLEQP_INACTIVE);

        SLEQP_CALL(sleqp_mat_push(augmented_matrix,
                                  num_variables + variable_index,
                                  column,
                                  1.));
      }
      else
      {
        assert(var_state == SLEQP_INACTIVE);
      }
    }

    int* cons_jac_cols    = sleqp_mat_cols(cons_jac);
    int* cons_jac_rows    = sleqp_mat_rows(cons_jac);
    double* cons_jac_data = sleqp_mat_data(cons_jac);

    for (int index = cons_jac_cols[column]; index < cons_jac_cols[column + 1];
         ++index)
    {
      int jac_row = cons_jac_rows[index];

      int act_cons_index = sleqp_working_set_cons_index(working_set, jac_row);

      SLEQP_ACTIVE_STATE cons_state
        = sleqp_working_set_cons_state(working_set, jac_row);

      if (act_cons_index != SLEQP_NONE)
      {
        assert(cons_state != SLEQP_INACTIVE);

        const double value = cons_jac_data[index];

        int aug_row = num_variables + act_cons_index;

        SLEQP_CALL(sleqp_mat_push(augmented_matrix,
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

  if (!lower_only)
  {
    SLEQP_CALL(add_upper(jacobian));
  }
  else
  {
    for (int j = num_variables; j < augmented_size; ++j)
    {
      SLEQP_CALL(sleqp_mat_push_col(augmented_matrix, j));
    }

    assert(sleqp_mat_is_lower(augmented_matrix));
  }
  // else: push columns!

  assert(sleqp_mat_is_valid(augmented_matrix));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
aug_jac_set_iterate(SleqpIterate* iterate, void* data)
{
  AugJacData* jacobian = (AugJacData*)data;

  SleqpProblem* problem        = jacobian->problem;
  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);

  // Do not recompute for linear problems & unchanged working set
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

  jacobian->condition = SLEQP_NONE;

  const bool lower_only
    = sleqp_fact_flags(jacobian->fact) & SLEQP_FACT_FLAGS_LOWER;

  SLEQP_CALL(reserve_aug_jac(jacobian, iterate, lower_only));
  SLEQP_CALL(fill_aug_jac(jacobian, iterate, lower_only));

  SLEQP_CALL(sleqp_timer_start(jacobian->factorization_timer));

  SLEQP_CALL(sleqp_fact_set_matrix(jacobian->fact, jacobian->augmented_matrix));

  jacobian->has_factorization = true;

  SLEQP_CALL(sleqp_timer_stop(jacobian->factorization_timer));

  SLEQP_CALL(sleqp_fact_cond(jacobian->fact, &jacobian->condition));

  {
    SleqpMat* matrix = jacobian->augmented_matrix;

    int total_size = num_variables + jacobian->working_set_size;

    assert(sleqp_mat_is_quadratic(matrix));
    assert(sleqp_mat_num_cols(matrix) == total_size);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
aug_jac_condition(bool* exact, double* condition, void* data)
{
  AugJacData* jacobian = (AugJacData*)data;

  *exact     = false;
  *condition = jacobian->condition;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
aug_jac_solve_min_norm(const SleqpVec* _rhs, SleqpVec* sol, void* data)
{
  AugJacData* jacobian = (AugJacData*)data;

  // Cast away constness
  SleqpVec* rhs = (SleqpVec*)_rhs;

  assert(jacobian->fact);

  SLEQP_CALL(sleqp_timer_start(jacobian->substitution_timer));

  SleqpProblem* problem    = jacobian->problem;
  SleqpFact* factorization = jacobian->fact;

  const double zero_eps
    = sleqp_params_value(jacobian->params, SLEQP_PARAM_ZERO_EPS);

  const int num_variables = sleqp_problem_num_vars(problem);

  assert(sol->dim == num_variables);

  rhs->dim += num_variables;

  for (int k = 0; k < rhs->nnz; ++k)
  {
    rhs->indices[k] += num_variables;
  }

  SLEQP_CALL(sleqp_fact_solve(factorization, rhs));

  SLEQP_CALL(
    sleqp_fact_solution(factorization, sol, 0, num_variables, zero_eps));

  rhs->dim -= num_variables;

  for (int k = 0; k < rhs->nnz; ++k)
  {
    rhs->indices[k] -= num_variables;
  }

  SLEQP_CALL(sleqp_timer_stop(jacobian->substitution_timer));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
aug_jac_solve_lsq(const SleqpVec* _rhs, SleqpVec* sol, void* data)
{
  AugJacData* jacobian = (AugJacData*)data;

  // Cast away constness
  SleqpVec* rhs = (SleqpVec*)_rhs;

  assert(jacobian->fact);

  SLEQP_CALL(sleqp_timer_start(jacobian->substitution_timer));

  SleqpProblem* problem    = jacobian->problem;
  SleqpFact* factorization = jacobian->fact;

  double zero_eps = sleqp_params_value(jacobian->params, SLEQP_PARAM_ZERO_EPS);

  const int num_variables    = sleqp_problem_num_vars(problem);
  const int working_set_size = jacobian->working_set_size;
  const int total_size       = num_variables + working_set_size;

  assert(rhs->dim == num_variables);

  // just add some zeros...
  SLEQP_CALL(sleqp_vec_resize(rhs, total_size));

  SLEQP_CALL(sleqp_fact_solve(factorization, rhs));

  assert(sol->dim == working_set_size);

  SLEQP_CALL(sleqp_fact_solution(factorization,
                                 sol,
                                 num_variables,
                                 total_size,
                                 zero_eps));

  // erase the zeros
  SLEQP_CALL(sleqp_vec_resize(rhs, num_variables));

  SLEQP_CALL(sleqp_timer_stop(jacobian->substitution_timer));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
aug_jac_project_nullspace(const SleqpVec* _rhs, SleqpVec* sol, void* data)
{
  AugJacData* jacobian = (AugJacData*)data;

  // Cast away constness
  SleqpVec* rhs = (SleqpVec*)_rhs;

  assert(jacobian->fact);

  SLEQP_CALL(sleqp_timer_start(jacobian->substitution_timer));

  SleqpProblem* problem    = jacobian->problem;
  SleqpFact* factorization = jacobian->fact;

  double zero_eps = sleqp_params_value(jacobian->params, SLEQP_PARAM_ZERO_EPS);

  const int num_variables    = sleqp_problem_num_vars(problem);
  const int working_set_size = jacobian->working_set_size;
  const int total_size       = num_variables + working_set_size;

  assert(rhs->dim == num_variables);

  // just add some zeros...
  SLEQP_CALL(sleqp_vec_resize(rhs, total_size));

  SLEQP_CALL(sleqp_fact_solve(factorization, rhs));

  assert(sol->dim == num_variables);

  SLEQP_CALL(
    sleqp_fact_solution(factorization, sol, 0, num_variables, zero_eps));

  // erase the zeros
  SLEQP_CALL(sleqp_vec_resize(rhs, num_variables));

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

  SLEQP_CALL(sleqp_fact_release(&jacobian->fact));

  SLEQP_CALL(sleqp_mat_release(&jacobian->augmented_matrix));

  SLEQP_CALL(sleqp_params_release(&jacobian->params));

  SLEQP_CALL(sleqp_problem_release(&jacobian->problem));

  sleqp_free(&jacobian);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
aug_jac_data_create(AugJacData** star,
                    SleqpProblem* problem,
                    SleqpParams* params,
                    SleqpFact* factorization)
{
  SLEQP_CALL(sleqp_malloc(star));

  const int num_variables   = sleqp_problem_num_vars(problem);
  const int num_constraints = sleqp_problem_num_cons(problem);

  AugJacData* jacobian = *star;

  *jacobian = (AugJacData){0};

  SLEQP_CALL(sleqp_problem_capture(problem));
  jacobian->problem = problem;

  SLEQP_CALL(sleqp_params_capture(params));
  jacobian->params = params;

  jacobian->condition = SLEQP_NONE;

  SLEQP_CALL(sleqp_mat_create(&jacobian->augmented_matrix, 0, 0, 0));

  jacobian->max_set_size = num_constraints + num_variables;
  int max_num_cols       = num_variables + jacobian->max_set_size;

  jacobian->has_factorization = false;

  SLEQP_CALL(sleqp_fact_capture(factorization));
  jacobian->fact = factorization;

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
                              SleqpFact* factorization)
{
  AugJacData* aug_jac_data;

  SLEQP_CALL(
    aug_jac_data_create(&aug_jac_data, problem, params, factorization));

  SleqpAugJacCallbacks callbacks
    = {.set_iterate       = aug_jac_set_iterate,
       .solve_min_norm    = aug_jac_solve_min_norm,
       .solve_lsq         = aug_jac_solve_lsq,
       .project_nullspace = aug_jac_project_nullspace,
       .condition         = aug_jac_condition,
       .free              = aug_jac_free};

  SLEQP_CALL(sleqp_aug_jac_create(star, problem, &callbacks, aug_jac_data));

  return SLEQP_OKAY;
}
