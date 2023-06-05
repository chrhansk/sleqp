#include "direct_aug_jac.h"

#include <assert.h>

#include "cmp.h"
#include "mem.h"

typedef struct
{
  SleqpFactQR* fact;

  SleqpVec* product;

  SleqpMat* working_jac;

  int* row_sums;
  // transposed working Jacobian
  SleqpMat* matrix;

} AugJacData;

static SLEQP_RETCODE
update_dims(AugJacData* aug_jac,
            const SleqpMat* cons_jac,
            SleqpWorkingSet* working_set)
{
  // store the transpose!
  const int cons_jac_num_rows = sleqp_mat_num_rows(cons_jac);
  const int cons_jac_num_cols = sleqp_mat_num_cols(cons_jac);

  SLEQP_CALL(sleqp_vec_clear(aug_jac->product));
  SLEQP_CALL(sleqp_vec_resize(aug_jac->product, cons_jac_num_rows));

  sleqp_free(&aug_jac->row_sums);
  SLEQP_CALL(sleqp_alloc_array(&aug_jac->row_sums, cons_jac_num_rows));

  const int num_active_cons = sleqp_working_set_num_active_cons(working_set);

  SLEQP_CALL(
    sleqp_mat_resize(aug_jac->working_jac, num_active_cons, cons_jac_num_cols));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_working_jac(AugJacData* aug_jac,
                    const SleqpMat* cons_jac,
                    SleqpWorkingSet* working_set)
{
  const int num_rows = sleqp_mat_num_rows(cons_jac);
  const int num_cols = sleqp_mat_num_cols(cons_jac);

  const int num_active_vars  = sleqp_working_set_num_active_vars(working_set);
  const int working_set_size = sleqp_working_set_size(working_set);

  int num_removed = 0;

  int* removed_rows = aug_jac->row_sums;

  for (int i = 0; i < num_rows; ++i)
  {
    if (sleqp_working_set_cons_state(working_set, i) == SLEQP_INACTIVE)
    {
      removed_rows[num_removed++] = i;
    }
  }

  SleqpMat* working_jac = aug_jac->working_jac;

  // Remove inactive cons
  SLEQP_CALL(
    sleqp_mat_remove_rows(cons_jac, working_jac, removed_rows, num_removed));

  if (num_active_vars == 0)
  {
    return SLEQP_OKAY;
  }

  int working_nnz = sleqp_mat_nnz(working_jac);
  working_nnz += num_active_vars;

  SLEQP_CALL(sleqp_mat_reserve(working_jac, working_nnz));

  SLEQP_CALL(sleqp_mat_resize(working_jac, working_set_size, num_cols));

  // Push active vars
  {
    int offset = num_active_vars;

    int* working_cols    = sleqp_mat_cols(working_jac);
    int* working_rows    = sleqp_mat_rows(working_jac);
    double* working_data = sleqp_mat_data(working_jac);

    for (int col = num_cols - 1; col >= 0; --col)
    {
      const int begin = working_cols[col];
      const int end   = working_cols[col + 1];

      for (int i = end - 1; i >= begin; --i)
      {
        working_rows[i + offset] = working_rows[i] + num_active_vars;
        working_data[i + offset] = working_data[i];
      }

      working_cols[col + 1] += offset;

      if (sleqp_working_set_var_state(working_set, col) != SLEQP_INACTIVE)
      {
        working_data[begin + offset - 1] = 1.;
        working_rows[begin + offset - 1] = offset - 1;

        if (--offset == 0)
        {
          break;
        }
      }
    }

    assert(offset == 0);
  }

  SLEQP_CALL(sleqp_mat_set_nnz(working_jac, working_nnz));

  assert(sleqp_mat_is_valid(working_jac));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
direct_aug_jac_set_iterate(SleqpIterate* iterate, void* data)
{
  AugJacData* aug_jac = (AugJacData*)data;

  SleqpMat* cons_jac           = sleqp_iterate_cons_jac(iterate);
  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);

  SLEQP_CALL(update_dims(aug_jac, cons_jac, working_set));

  SLEQP_CALL(compute_working_jac(aug_jac, cons_jac, working_set));

  SleqpMat* working_jac = aug_jac->working_jac;

  int* row_sums = aug_jac->row_sums;

  SLEQP_CALL(sleqp_mat_trans(working_jac, aug_jac->matrix, row_sums));

  assert(sleqp_mat_is_valid(aug_jac->matrix));

  SLEQP_CALL(sleqp_qr_set_matrix(aug_jac->fact, aug_jac->matrix));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
direct_aug_jac_solve_min_norm(const SleqpVec* rhs, SleqpVec* sol, void* data)
{
  AugJacData* aug_jac = (AugJacData*)data;

  SleqpMat* matrix  = aug_jac->matrix;
  SleqpVec* product = aug_jac->product;

  const int num_rows = sleqp_mat_num_rows(matrix);
  const int num_cols = sleqp_mat_num_cols(matrix);

  SLEQP_CALL(sleqp_vec_clear(product));
  SLEQP_CALL(sleqp_vec_resize(product, num_cols));

  SLEQP_CALL(sleqp_qr_solve_tri_trans(aug_jac->fact, rhs, product));

  // Enlarge by adding zeros for the null space part
  SLEQP_CALL(sleqp_vec_resize(product, num_rows));

  SLEQP_CALL(sleqp_qr_mult_orth(aug_jac->fact, product, sol));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
direct_aug_jac_solve_lsq(const SleqpVec* rhs, SleqpVec* sol, void* data)
{
  AugJacData* aug_jac = (AugJacData*)data;

  SleqpMat* matrix  = aug_jac->matrix;
  SleqpVec* product = aug_jac->product;

  const int num_rows = sleqp_mat_num_rows(matrix);
  const int num_cols = sleqp_mat_num_cols(matrix);

  assert(rhs->dim == num_rows);

  SLEQP_CALL(sleqp_vec_resize(product, num_rows));

  SLEQP_CALL(sleqp_qr_mult_orth_trans(aug_jac->fact, rhs, product));

  // Discard elements associated with null space part
  SLEQP_CALL(sleqp_vec_resize(product, num_cols));

  SLEQP_CALL(sleqp_qr_solve_tri(aug_jac->fact, product, sol));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
direct_aug_jac_project_nullspace(const SleqpVec* rhs, SleqpVec* sol, void* data)
{
  AugJacData* aug_jac = (AugJacData*)data;

  SleqpMat* matrix  = aug_jac->matrix;
  SleqpVec* product = aug_jac->product;

  const int num_rows = sleqp_mat_num_rows(matrix);
  const int num_cols = sleqp_mat_num_cols(matrix);

  product->dim = num_rows;

  SLEQP_CALL(sleqp_qr_mult_orth_trans(aug_jac->fact, rhs, product));

  // slice off first num_cols entries
  {
    int k = 0;

    for (; k < product->nnz; ++k)
    {
      if (product->indices[k] >= num_cols)
      {
        break;
      }
    }

    const int offset = k;

    for (; k < product->nnz; ++k)
    {
      product->indices[k - offset] = product->indices[k];
      product->data[k - offset]    = product->data[k];
    }

    product->nnz -= offset;
  }

  SLEQP_CALL(sleqp_qr_mult_orth(aug_jac->fact, product, sol));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
direct_aug_jac_condition(bool* exact, double* condition, void* data)
{
  *condition = SLEQP_NONE;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
direct_aug_jac_free(void* data)
{
  AugJacData* aug_jac = data;

  SLEQP_CALL(sleqp_qr_release(&aug_jac->fact));

  SLEQP_CALL(sleqp_mat_release(&aug_jac->matrix));

  sleqp_free(&aug_jac->row_sums);

  SLEQP_CALL(sleqp_mat_release(&aug_jac->working_jac));

  SLEQP_CALL(sleqp_vec_free(&aug_jac->product));

  sleqp_free(&data);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
create_aug_jac_data(AugJacData** star, SleqpFactQR* fact)
{
  SLEQP_CALL(sleqp_malloc(star));

  AugJacData* aug_jac = *star;

  *aug_jac = (AugJacData){};

  SLEQP_CALL(sleqp_vec_create_empty(&aug_jac->product, 0));

  SLEQP_CALL(sleqp_mat_create(&aug_jac->working_jac, 0, 0, 0));

  SLEQP_CALL(sleqp_mat_create(&aug_jac->matrix, 0, 0, 0));
  sleqp_free(&aug_jac->row_sums);

  SLEQP_CALL(sleqp_qr_capture(fact));
  aug_jac->fact = fact;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_direct_aug_jac_create(SleqpAugJac** star,
                            SleqpProblem* problem,
                            SleqpSettings* settings,
                            SleqpFactQR* fact)
{
  AugJacData* data = NULL;

  SLEQP_CALL(create_aug_jac_data(&data, fact));

  SleqpAugJacCallbacks callbacks = {
    .set_iterate       = direct_aug_jac_set_iterate,
    .solve_min_norm    = direct_aug_jac_solve_min_norm,
    .solve_lsq         = direct_aug_jac_solve_lsq,
    .project_nullspace = direct_aug_jac_project_nullspace,
    .condition         = direct_aug_jac_condition,
    .free              = direct_aug_jac_free,
  };

  SLEQP_CALL(sleqp_aug_jac_create(star, problem, &callbacks, data));

  return SLEQP_OKAY;
}
