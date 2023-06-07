#include "reduced_aug_jac.h"

#include "aug_jac/aug_jac.h"
#include "cmp.h"
#include "fact/fact.h"
#include "iterate.h"
#include "mem.h"
#include "problem.h"
#include "sparse/mat.h"
#include "working_set.h"

typedef struct
{
  SleqpProblem* problem;
  SleqpSettings* settings;

  SleqpFact* fact;
  bool has_factorization;

  // J^T*J
  SleqpMat* matrix;

  SleqpVec* rhs;
  SleqpVec* direction;

  int num_rows;
  int num_working_cons;
  int num_nnz;

  int* row_counts;
  int* row_offsets;

  int* col_indices;
  double* values;

  double* cache;

  SleqpWorkingSet* working_set;

} AugJacData;

static SLEQP_RETCODE
aug_jac_data_create(AugJacData** star,
                    SleqpProblem* problem,
                    SleqpSettings* settings,
                    SleqpFact* factorization)
{
  SLEQP_CALL(sleqp_malloc(star));

  AugJacData* jacobian = *star;

  *jacobian = (AugJacData){0};

  SLEQP_CALL(sleqp_problem_capture(problem));
  jacobian->problem = problem;

  SLEQP_CALL(sleqp_settings_capture(settings));
  jacobian->settings = settings;

  SLEQP_CALL(sleqp_fact_capture(factorization));
  jacobian->fact              = factorization;
  jacobian->has_factorization = false;

  SLEQP_CALL(sleqp_mat_create(&jacobian->matrix, 0, 0, 0));

  SLEQP_CALL(sleqp_vec_create_empty(&jacobian->rhs, 0));

  const int num_vars = sleqp_problem_num_vars(problem);

  SLEQP_CALL(sleqp_vec_create_full(&jacobian->direction, num_vars));

  SLEQP_CALL(sleqp_alloc_array(&jacobian->cache, num_vars));

  const bool fixed_jacobian = !(sleqp_problem_has_nonlinear_cons(problem));

  if (fixed_jacobian)
  {
    SLEQP_CALL(sleqp_working_set_create(&jacobian->working_set, problem));
  }

  jacobian->num_rows         = SLEQP_NONE;
  jacobian->num_working_cons = SLEQP_NONE;
  jacobian->num_nnz          = SLEQP_NONE;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
reserve_rows(AugJacData* jacobian, SleqpWorkingSet* working_set)
{
  const int num_rows = sleqp_working_set_size(working_set);

  if (jacobian->num_rows != num_rows)
  {
    SLEQP_CALL(sleqp_realloc(&jacobian->row_counts, num_rows + 1));
    SLEQP_CALL(sleqp_vec_resize(jacobian->rhs, num_rows));
    SLEQP_CALL(sleqp_vec_reserve(jacobian->rhs, num_rows));
    jacobian->num_rows = num_rows;
  }

  const int num_working_cons = sleqp_working_set_num_active_cons(working_set);

  if (jacobian->num_working_cons != num_working_cons)
  {
    SLEQP_CALL(sleqp_realloc(&jacobian->row_offsets, num_working_cons));
    jacobian->num_working_cons = num_working_cons;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_row_counts_and_nnz(AugJacData* jacobian,
                           SleqpMat* cons_jac,
                           SleqpWorkingSet* working_set,
                           int* nnz)
{
  *nnz = sleqp_working_set_num_active_vars(working_set);

  const int num_active_vars = sleqp_working_set_num_active_vars(working_set);

  for (int i = 0; i < num_active_vars; ++i)
  {
    jacobian->row_counts[i] = 1;
  }

  assert(sleqp_working_set_size(working_set) == jacobian->num_rows);

  for (int i = num_active_vars; i < jacobian->num_rows; ++i)
  {
    jacobian->row_counts[i] = 0;
  }

  const int num_cols = sleqp_mat_num_cols(cons_jac);

  const int* cols = sleqp_mat_cols(cons_jac);
  const int* rows = sleqp_mat_rows(cons_jac);

  for (int col = 0; col < num_cols; ++col)
  {
    for (int index = cols[col]; index < cols[col + 1]; ++index)
    {
      const int row        = rows[index];
      const int cons_index = sleqp_working_set_cons_index(working_set, row);

      if (cons_index == SLEQP_NONE)
      {
        continue;
      }

      assert(cons_index >= num_active_vars);
      assert(cons_index < jacobian->num_rows);

      ++(jacobian->row_counts[cons_index]);
      ++(*nnz);
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
reserve_data(AugJacData* jacobian, int num_nnz)
{
  if (jacobian->num_nnz != num_nnz)
  {
    SLEQP_CALL(sleqp_realloc(&jacobian->col_indices, num_nnz));
    SLEQP_CALL(sleqp_realloc(&jacobian->values, num_nnz));
    jacobian->num_nnz = num_nnz;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
update_row_counts(AugJacData* jacobian)
{
  int offset = 0;

  for (int i = 0; i <= jacobian->num_rows; ++i)
  {
    const int current       = jacobian->row_counts[i];
    jacobian->row_counts[i] = offset;
    offset += current;
  }

  for (int i = 0; i < jacobian->num_working_cons; ++i)
  {
    jacobian->row_offsets[i] = 0;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
fill_system(AugJacData* jacobian,
            SleqpMat* cons_jac,
            SleqpWorkingSet* working_set)
{
  const int num_active_vars = sleqp_working_set_num_active_vars(working_set);

  {
    for (int i = 0; i < num_active_vars; ++i)
    {
      jacobian->values[i] = 1.;
    }

    SleqpProblem* problem = jacobian->problem;

    const int num_vars = sleqp_problem_num_vars(problem);

    for (int j = 0; j < num_vars; ++j)
    {
      const int index = sleqp_working_set_var_index(working_set, j);

      if (index == SLEQP_NONE)
      {
        continue;
      }

      jacobian->col_indices[index] = j;
    }
  }

  {
    const int num_cols   = sleqp_mat_num_cols(cons_jac);
    const int* cols      = sleqp_mat_cols(cons_jac);
    const int* rows      = sleqp_mat_rows(cons_jac);
    const double* values = sleqp_mat_data(cons_jac);

    for (int col = 0; col < num_cols; ++col)
    {
      for (int index = cols[col]; index < cols[col + 1]; ++index)
      {
        const int row  = rows[index];
        int cons_index = sleqp_working_set_cons_index(working_set, row);

        if (cons_index == SLEQP_NONE)
        {
          continue;
        }

        assert(cons_index >= num_active_vars);

        cons_index -= num_active_vars;

        const int jac_index = jacobian->row_counts[cons_index + num_active_vars]
                              + jacobian->row_offsets[cons_index];

        jacobian->values[jac_index]      = values[index];
        jacobian->col_indices[jac_index] = col;

        ++(jacobian->row_offsets[cons_index]);
      }
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
reserve_more_and_push(SleqpMat* matrix, int row, int col, double value)
{
  const int nnz = sleqp_mat_nnz(matrix);

  if (nnz == sleqp_mat_nnz_max(matrix))
  {
    const int num_cols = sleqp_mat_num_cols(matrix);
    const int num_rows = sleqp_mat_num_rows(matrix);

    const int size = num_cols * num_rows;

    const int next_nnz = SLEQP_MIN(size, 2 * nnz + 1);
    SLEQP_CALL(sleqp_mat_reserve(matrix, next_nnz));
  }

  SLEQP_CALL(sleqp_mat_push(matrix, row, col, value));

  return SLEQP_OKAY;
}

void
compute_inner_product(const int* rows,
                      const int* cols,
                      const double* data,
                      int row,
                      int col,
                      bool* nonzero,
                      double* product)
{
  *nonzero = true;
  *product = 0.;

  int col_index = rows[col];
  int row_index = rows[row];

  const int col_bound = rows[col + 1];
  const int row_bound = rows[row + 1];

  while (col_index < col_bound && row_index < row_bound)
  {
    const int col_col = cols[col_index];
    const int row_col = cols[row_index];

    if (row_col < col_col)
    {
      ++row_index;
    }
    else if (row_col > col_col)
    {
      ++col_index;
    }
    else
    {
      (*nonzero) = true;
      (*product) += data[row_index] * data[col_index];
      ++row_index;
      ++col_index;
    }
  }
}

static SLEQP_RETCODE
compute_matrix_lower(AugJacData* jacobian, int num_active_vars)
{
  const int size = jacobian->num_rows;

  const int* rows    = jacobian->row_counts;
  const int* cols    = jacobian->col_indices;
  const double* data = jacobian->values;

  SleqpMat* matrix = jacobian->matrix;

  SLEQP_CALL(sleqp_mat_clear(matrix));
  SLEQP_CALL(sleqp_mat_resize(matrix, size, size));

  for (int col = 0; col < size; ++col)
  {
    SLEQP_CALL(sleqp_mat_push_col(matrix, col));

    // diagonal

    if (col < num_active_vars)
    {
      SLEQP_CALL(reserve_more_and_push(matrix, col, col, 1.));
    }
    else
    {
      double product = 0.;

      for (int col_index = rows[col]; col_index < rows[col + 1]; ++col_index)
      {
        product += data[col_index] * data[col_index];
      }

      SLEQP_CALL(reserve_more_and_push(matrix, col, col, product));
    }

    for (int row = col + 1; row < size; ++row)
    {
      double product;
      bool nonzero;

      compute_inner_product(rows, cols, data, row, col, &nonzero, &product);

      // push product
      if (nonzero)
      {
        SLEQP_CALL(reserve_more_and_push(matrix, row, col, product));
      }
    }
  }

  assert(sleqp_mat_is_lower(matrix));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_matrix_general(AugJacData* jacobian)
{
  const int size = jacobian->num_rows;

  const int* rows    = jacobian->row_counts;
  const int* cols    = jacobian->col_indices;
  const double* data = jacobian->values;

  SleqpMat* matrix = jacobian->matrix;

  SLEQP_CALL(sleqp_mat_clear(matrix));
  SLEQP_CALL(sleqp_mat_resize(matrix, size, size));

  for (int col = 0; col < size; ++col)
  {
    SLEQP_CALL(sleqp_mat_push_col(matrix, col));

    for (int row = 0; row < size; ++row)
    {
      double product;
      bool nonzero;

      compute_inner_product(rows, cols, data, row, col, &nonzero, &product);

      // push product
      if (nonzero)
      {
        SLEQP_CALL(reserve_more_and_push(matrix, row, col, product));
      }
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_matrix(AugJacData* jacobian, int num_active_vars)
{
  SleqpFact* factorization = jacobian->fact;

  if (sleqp_fact_flags(factorization) & SLEQP_FACT_FLAGS_LOWER)
  {
    SLEQP_CALL(compute_matrix_lower(jacobian, num_active_vars));
  }
  else
  {
    SLEQP_CALL(compute_matrix_general(jacobian));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
aug_jac_set_iterate(SleqpIterate* iterate, void* data)
{
  AugJacData* jacobian = (AugJacData*)data;

  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);
  SleqpMat* cons_jac           = sleqp_iterate_cons_jac(iterate);

  const int num_active_vars = sleqp_working_set_num_active_vars(working_set);

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

  SLEQP_CALL(reserve_rows(jacobian, working_set));

  int total_nnz = SLEQP_NONE;

  SLEQP_CALL(compute_row_counts_and_nnz(jacobian,
                                        sleqp_iterate_cons_jac(iterate),
                                        working_set,
                                        &total_nnz));

  SLEQP_CALL(reserve_data(jacobian, total_nnz));

  SLEQP_CALL(update_row_counts(jacobian));

  SLEQP_CALL(fill_system(jacobian, cons_jac, working_set));

  SLEQP_CALL(compute_matrix(jacobian, num_active_vars));

  SLEQP_CALL(sleqp_fact_set_matrix(jacobian->fact, jacobian->matrix));
  jacobian->has_factorization = true;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_product(AugJacData* jacobian,
                const SleqpVec* direction,
                SleqpVec* product)
{
  const int size = jacobian->num_rows;

  assert(direction != product);

  SLEQP_CALL(sleqp_vec_clear(product));
  SLEQP_CALL(sleqp_vec_reserve(product, size));

  const int* rows = jacobian->row_counts;

  for (int row = 0; row < size; ++row)
  {
    int mat_index       = rows[row];
    const int mat_bound = rows[row + 1];

    double prod  = 0.;
    bool nonzero = false;

    int rhs_index = 0;

    while (mat_index < mat_bound && rhs_index < direction->nnz)
    {
      const int mat_col = jacobian->col_indices[mat_index];
      const int vec_col = direction->indices[rhs_index];

      if (mat_col < vec_col)
      {
        ++mat_index;
      }
      else if (mat_col > vec_col)
      {
        ++rhs_index;
      }
      else
      {
        nonzero = true;
        prod += jacobian->values[mat_index] * direction->data[rhs_index];

        ++mat_index;
        ++rhs_index;
      }
    }

    if (nonzero)
    {
      SLEQP_CALL(sleqp_vec_push(product, row, prod));
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_trans_product(AugJacData* jacobian, const SleqpVec* direction)
{
  SleqpProblem* problem = jacobian->problem;

  const int num_vars = sleqp_problem_num_vars(problem);

  for (int i = 0; i < num_vars; ++i)
  {
    jacobian->cache[i] = 0.;
  }

  const int* rows = jacobian->row_counts;
  const int* cols = jacobian->col_indices;

  for (int k = 0; k < direction->nnz; ++k)
  {
    const int row      = direction->indices[k];
    const double value = direction->data[k];

    for (int index = rows[row]; index < rows[row + 1]; ++index)
    {
      int col = cols[index];
      jacobian->cache[col] += value * jacobian->values[index];
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
aug_jac_solve_min_norm(const SleqpVec* rhs, SleqpVec* sol, void* data)
{
  AugJacData* jacobian  = (AugJacData*)data;
  SleqpProblem* problem = jacobian->problem;

  const int num_vars = sleqp_problem_num_vars(problem);
  const int size     = jacobian->num_rows;

  SleqpVec* product = jacobian->rhs;

  const double zero_eps
    = sleqp_settings_real_value(jacobian->settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  SLEQP_CALL(sleqp_fact_solve(jacobian->fact, rhs));

  SLEQP_CALL(sleqp_fact_solution(jacobian->fact, product, 0, size, zero_eps));

  SLEQP_CALL(compute_trans_product(jacobian, product));

  SLEQP_CALL(sleqp_vec_set_from_raw(sol, jacobian->cache, num_vars, zero_eps));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
aug_jac_solve_lsq(const SleqpVec* rhs, SleqpVec* sol, void* data)
{
  AugJacData* jacobian = (AugJacData*)data;

  const int size = jacobian->num_rows;

  SleqpVec* product = jacobian->rhs;

  const double zero_eps
    = sleqp_settings_real_value(jacobian->settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  SLEQP_CALL(compute_product(jacobian, rhs, product));

  SLEQP_CALL(sleqp_fact_solve(jacobian->fact, product));

  SLEQP_CALL(sleqp_fact_solution(jacobian->fact, sol, 0, size, zero_eps));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
aug_jac_project_nullspace(const SleqpVec* rhs, SleqpVec* sol, void* data)
{
  AugJacData* jacobian  = (AugJacData*)data;
  SleqpProblem* problem = jacobian->problem;

  const int num_vars = sleqp_problem_num_vars(problem);

  const int size = jacobian->num_rows;

  SleqpVec* product = jacobian->rhs;

  const double zero_eps
    = sleqp_settings_real_value(jacobian->settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  SLEQP_CALL(compute_product(jacobian, rhs, product));

  SLEQP_CALL(sleqp_fact_solve(jacobian->fact, product));

  SLEQP_CALL(sleqp_fact_solution(jacobian->fact, product, 0, size, zero_eps));

  SLEQP_CALL(compute_trans_product(jacobian, product));

  SLEQP_CALL(sleqp_vec_set_from_raw(jacobian->direction,
                                    jacobian->cache,
                                    num_vars,
                                    zero_eps));

  SLEQP_CALL(
    sleqp_vec_add_scaled(rhs, jacobian->direction, 1., -1., zero_eps, sol));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
aug_jac_condition(bool* exact, double* condition, void* data)
{
  assert(false);
  return SLEQP_OKAY;
}

static SLEQP_RETCODE
aug_jac_free(void* data)
{
  AugJacData* jacobian = (AugJacData*)data;

  SLEQP_CALL(sleqp_working_set_release(&jacobian->working_set));

  sleqp_free(&jacobian->cache);
  sleqp_free(&jacobian->values);

  sleqp_free(&jacobian->col_indices);
  sleqp_free(&jacobian->row_offsets);
  sleqp_free(&jacobian->row_counts);

  SLEQP_CALL(sleqp_vec_free(&jacobian->direction));
  SLEQP_CALL(sleqp_vec_free(&jacobian->rhs));
  SLEQP_CALL(sleqp_mat_release(&jacobian->matrix));

  SLEQP_CALL(sleqp_fact_release(&jacobian->fact));

  SLEQP_CALL(sleqp_settings_release(&jacobian->settings));

  SLEQP_CALL(sleqp_problem_release(&jacobian->problem));

  sleqp_free(&jacobian);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_reduced_aug_jac_create(SleqpAugJac** star,
                             SleqpProblem* problem,
                             SleqpSettings* settings,
                             SleqpFact* factorization)
{
  AugJacData* aug_jac_data;

  SLEQP_CALL(
    aug_jac_data_create(&aug_jac_data, problem, settings, factorization));

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
