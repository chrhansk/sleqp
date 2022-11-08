#include "fixed_var_func.h"

#include "cmp.h"
#include "dyn.h"
#include "fail.h"
#include "lsq.h"
#include "mem.h"

#include "preprocessor/preprocessing.h"

typedef struct FixedVarFuncData
{
  int num_fixed;
  double* fixed_values;
  int* fixed_indices;

  SleqpFunc* func;

  SleqpVec* values;

  SleqpVec* grad;

  SleqpVec* direction;
  SleqpVec* product;

  SleqpMat* jacobian;

} FixedVarFuncData;

static SLEQP_RETCODE
fixed_var_func_set(SleqpFunc* func,
                   SleqpVec* value,
                   SLEQP_VALUE_REASON reason,
                   bool* reject,
                   void* data)
{
  FixedVarFuncData* func_data = (FixedVarFuncData*)data;

  SLEQP_CALL(sleqp_preprocessing_merge_entries(value,
                                               func_data->values,
                                               func_data->num_fixed,
                                               func_data->fixed_indices,
                                               func_data->fixed_values));

  SLEQP_CALL(
    sleqp_func_set_value(func_data->func, func_data->values, reason, reject));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
fixed_var_func_nonzeros(SleqpFunc* func,
                        int* obj_grad_nnz,
                        int* cons_val_nnz,
                        int* cons_jac_nnz,
                        int* hess_prod_nnz,
                        void* data)
{
  const int num_vars = sleqp_func_num_vars(func);
  const int num_cons = sleqp_func_num_cons(func);

  FixedVarFuncData* func_data = (FixedVarFuncData*)data;

  const int num_orig_vars = sleqp_func_num_vars(func_data->func);

  SLEQP_CALL(sleqp_func_nonzeros(func_data->func,
                                 obj_grad_nnz,
                                 cons_val_nnz,
                                 cons_jac_nnz,
                                 hess_prod_nnz));

  if (*obj_grad_nnz != SLEQP_NONE)
  {
    SLEQP_CALL(sleqp_vec_reserve(func_data->grad, *obj_grad_nnz));
  }
  else
  {
    SLEQP_CALL(sleqp_vec_reserve(func_data->grad, num_orig_vars));
  }

  if (*cons_jac_nnz != SLEQP_NONE)
  {
    SLEQP_CALL(sleqp_mat_reserve(func_data->jacobian, *cons_jac_nnz));
  }
  else
  {
    const int jac_nnz = num_orig_vars * num_cons;
    SLEQP_CALL(sleqp_mat_reserve(func_data->jacobian, jac_nnz));
  }

  if (*obj_grad_nnz != SLEQP_NONE)
  {
    *obj_grad_nnz = SLEQP_MIN(num_vars, *obj_grad_nnz);
  }

  if (*hess_prod_nnz != SLEQP_NONE)
  {
    *hess_prod_nnz = SLEQP_MIN(num_vars, *hess_prod_nnz);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
fixed_var_obj_val(SleqpFunc* func, double* obj_val, void* data)
{
  FixedVarFuncData* func_data = (FixedVarFuncData*)data;

  return sleqp_func_obj_val(func_data->func, obj_val);
}

static SLEQP_RETCODE
fixed_var_obj_grad(SleqpFunc* func, SleqpVec* obj_grad, void* data)
{
  FixedVarFuncData* func_data = (FixedVarFuncData*)data;

  SLEQP_CALL(sleqp_func_obj_grad(func_data->func, func_data->grad));

  SLEQP_CALL(sleqp_vec_remove_entries(func_data->grad,
                                      obj_grad,
                                      func_data->fixed_indices,
                                      func_data->num_fixed));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
fixed_var_cons_val(SleqpFunc* func, SleqpVec* cons_val, void* data)
{
  FixedVarFuncData* func_data = (FixedVarFuncData*)data;

  return sleqp_func_cons_val(func_data->func, cons_val);
}

static SLEQP_RETCODE
fixed_var_cons_jac(SleqpFunc* func, SleqpMat* cons_jac, void* data)
{
  FixedVarFuncData* func_data = (FixedVarFuncData*)data;

  SLEQP_CALL(sleqp_func_cons_jac(func_data->func, func_data->jacobian));

  SLEQP_CALL(sleqp_mat_remove_cols(func_data->jacobian,
                                   cons_jac,
                                   func_data->fixed_indices,
                                   func_data->num_fixed));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
fixed_var_hess_prod(SleqpFunc* func,
                    const double* obj_dual,
                    const SleqpVec* direction,
                    const SleqpVec* cons_duals,
                    SleqpVec* product,
                    void* data)
{
  FixedVarFuncData* func_data = (FixedVarFuncData*)data;

  SLEQP_CALL(sleqp_preprocessing_add_zero_entries(direction,
                                                  func_data->direction,
                                                  func_data->num_fixed,
                                                  func_data->fixed_indices));

  SLEQP_CALL(sleqp_func_hess_prod(func_data->func,
                                  obj_dual,
                                  func_data->direction,
                                  cons_duals,
                                  func_data->product));

  SLEQP_CALL(sleqp_vec_remove_entries(func_data->product,
                                      product,
                                      func_data->fixed_indices,
                                      func_data->num_fixed));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
fixed_func_free(void* data)
{
  FixedVarFuncData* func_data = (FixedVarFuncData*)data;

  SLEQP_CALL(sleqp_mat_release(&func_data->jacobian));

  SLEQP_CALL(sleqp_vec_free(&func_data->product));

  SLEQP_CALL(sleqp_vec_free(&func_data->direction));

  SLEQP_CALL(sleqp_vec_free(&func_data->grad));

  SLEQP_CALL(sleqp_vec_free(&func_data->values));

  SLEQP_CALL(sleqp_func_release(&func_data->func));

  sleqp_free(&func_data->fixed_values);

  sleqp_free(&func_data->fixed_indices);

  sleqp_free(&func_data);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
create_fixed_var_hess_struct(const SleqpHessStruct* source,
                             SleqpHessStruct* target,
                             const int num_fixed,
                             const int* fixed_indices)
{
  SLEQP_CALL(sleqp_hess_struct_clear(target));

  const int num_blocks = sleqp_hess_struct_num_blocks(source);

  int fixed_pos = 0;

  int offset = 0;
  for (int block = 0; block < num_blocks; ++block)
  {
    int source_begin, source_end;

    int next_offset = offset;

    SLEQP_CALL(
      sleqp_hess_struct_block_range(source, block, &source_begin, &source_end));

    const int target_begin = source_begin - offset;

    for (int j = source_begin; j < source_end; ++j)
    {
      if (fixed_pos < num_fixed && fixed_indices[fixed_pos] == j)
      {
        ++next_offset;
        ++fixed_pos;
      }
    }

    const int target_end = source_end - next_offset;

    if (target_begin != target_end)
    {
      SLEQP_CALL(sleqp_hess_struct_push_block(target, target_end));
    }

    offset = next_offset;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
fixed_lsq_func_residuals(SleqpFunc* func, SleqpVec* residuals, void* data)
{
  FixedVarFuncData* func_data = (FixedVarFuncData*)data;

  return sleqp_lsq_func_residuals(func_data->func, residuals);
}

static SLEQP_RETCODE
fixed_lsq_func_jac_forward(SleqpFunc* func,
                           const SleqpVec* forward_direction,
                           SleqpVec* product,
                           void* data)
{
  FixedVarFuncData* func_data = (FixedVarFuncData*)data;

  SLEQP_CALL(sleqp_preprocessing_add_zero_entries(forward_direction,
                                                  func_data->direction,
                                                  func_data->num_fixed,
                                                  func_data->fixed_indices));

  SLEQP_CALL(
    sleqp_lsq_func_jac_forward(func_data->func, func_data->direction, product));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
fixed_lsq_func_nonzeros(SleqpFunc* func,
                        int* residual_nnz,
                        int* jac_fwd_nnz,
                        int* jac_adj_nnz,
                        int* cons_val_nnz,
                        int* cons_jac_nnz,
                        void* data)
{
  FixedVarFuncData* func_data = (FixedVarFuncData*)data;

  const int num_variables = sleqp_func_num_vars(func);

  SLEQP_CALL(sleqp_lsq_func_nonzeros(func_data->func,
                                     residual_nnz,
                                     jac_fwd_nnz,
                                     jac_adj_nnz,
                                     cons_val_nnz,
                                     cons_jac_nnz));

  if (*cons_jac_nnz != SLEQP_NONE)
  {
    SLEQP_CALL(sleqp_mat_reserve(func_data->jacobian, *cons_jac_nnz));
  }

  if (*jac_adj_nnz != SLEQP_NONE)
  {
    *jac_adj_nnz = SLEQP_MIN(num_variables, *jac_adj_nnz);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
fixed_lsq_func_jac_adjoint(SleqpFunc* func,
                           const SleqpVec* adjoint_direction,
                           SleqpVec* product,
                           void* data)
{
  FixedVarFuncData* func_data = (FixedVarFuncData*)data;

  SleqpVec* full_product = func_data->direction;

  SLEQP_CALL(sleqp_lsq_func_jac_adjoint(func_data->func,
                                        adjoint_direction,
                                        full_product));

  SLEQP_CALL(sleqp_vec_remove_entries(full_product,
                                      product,
                                      func_data->fixed_indices,
                                      func_data->num_fixed));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
fixed_dyn_func_eval(SleqpFunc* func,
                    double* obj_val,
                    SleqpVec* cons_val,
                    double* error,
                    void* data)
{
  FixedVarFuncData* func_data = (FixedVarFuncData*)data;

  SLEQP_CALL(sleqp_dyn_func_eval(func_data->func, obj_val, cons_val, error));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
fixed_dyn_func_set_error_bound(SleqpFunc* func, double error_bound, void* data)
{
  FixedVarFuncData* func_data = (FixedVarFuncData*)data;

  SLEQP_CALL(sleqp_dyn_func_set_error_bound(func_data->func, error_bound));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
fixed_dyn_func_set_obj_weight(SleqpFunc* func, double obj_weight, void* data)
{
  FixedVarFuncData* func_data = (FixedVarFuncData*)data;

  SLEQP_CALL(sleqp_dyn_func_set_obj_weight(func_data->func, obj_weight));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
fixed_dyn_func_set_cons_weights(SleqpFunc* func,
                                const double* cons_weights,
                                void* data)
{
  FixedVarFuncData* func_data = (FixedVarFuncData*)data;

  SLEQP_CALL(sleqp_dyn_func_set_cons_weights(func_data->func, cons_weights));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
create_fixed_var_func_data(FixedVarFuncData** star,
                           SleqpFunc* func,
                           int num_fixed,
                           const int* fixed_indices,
                           const double* fixed_values)
{
  const int num_variables   = sleqp_func_num_vars(func);
  const int num_constraints = sleqp_func_num_cons(func);

  assert(num_fixed >= 0);
  assert(num_fixed <= num_variables);

  SLEQP_CALL(sleqp_malloc(star));

  FixedVarFuncData* func_data = *star;

  *func_data = (FixedVarFuncData){0};

  func_data->num_fixed = num_fixed;

  SLEQP_CALL(sleqp_alloc_array(&func_data->fixed_indices, num_fixed));

  SLEQP_CALL(sleqp_alloc_array(&func_data->fixed_values, num_fixed));

  int last_index = -1;

  for (int j = 0; j < num_fixed; ++j)
  {
    const int fixed_index = fixed_indices[j];

    sleqp_assert(fixed_index >= 0);
    sleqp_assert(fixed_index < num_variables);
    sleqp_assert(fixed_index > last_index);

    func_data->fixed_indices[j] = fixed_index;

    last_index = fixed_index;
  }

  for (int j = 0; j < num_fixed; ++j)
  {
    const double fixed_value = fixed_values[j];

    sleqp_assert(sleqp_is_finite(fixed_value));

    func_data->fixed_values[j] = fixed_value;
  }

  func_data->func = func;

  SLEQP_CALL(sleqp_func_capture(func));

  SLEQP_CALL(sleqp_vec_create_full(&func_data->values, num_variables));

  SLEQP_CALL(sleqp_vec_create_full(&func_data->grad, num_variables));

  SLEQP_CALL(sleqp_vec_create_full(&func_data->direction, num_variables));

  SLEQP_CALL(sleqp_vec_create_full(&func_data->product, num_variables));

  SLEQP_CALL(
    sleqp_mat_create(&func_data->jacobian, num_constraints, num_variables, 0));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_fixed_var_func_create(SleqpFunc** star,
                            SleqpFunc* func,
                            int num_fixed,
                            const int* fixed_indices,
                            const double* fixed_values)
{
  FixedVarFuncData* func_data;

  assert(sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_REGULAR);

  const int num_variables   = sleqp_func_num_vars(func);
  const int num_constraints = sleqp_func_num_cons(func);

  assert(num_fixed >= 0);
  assert(num_fixed <= num_variables);

  SLEQP_CALL(create_fixed_var_func_data(&func_data,
                                        func,
                                        num_fixed,
                                        fixed_indices,
                                        fixed_values));

  SleqpFuncCallbacks callbacks = {.set_value = fixed_var_func_set,
                                  .nonzeros  = fixed_var_func_nonzeros,
                                  .obj_val   = fixed_var_obj_val,
                                  .obj_grad  = fixed_var_obj_grad,
                                  .cons_val  = fixed_var_cons_val,
                                  .cons_jac  = fixed_var_cons_jac,
                                  .hess_prod = fixed_var_hess_prod,
                                  .func_free = fixed_func_free};

  SLEQP_CALL(sleqp_func_create(star,
                               &callbacks,
                               num_variables - num_fixed,
                               num_constraints,
                               (void*)func_data));

  SleqpFunc* fixed_var_func = *star;

  SLEQP_CALL(
    sleqp_func_flags_copy(func,
                          fixed_var_func,
                          SLEQP_FUNC_HESS_INEXACT | SLEQP_FUNC_HESS_PSD));

  SLEQP_CALL(sleqp_func_flags_add(fixed_var_func, SLEQP_FUNC_INTERNAL));

  SLEQP_CALL(
    create_fixed_var_hess_struct(sleqp_func_hess_struct(func),
                                 sleqp_func_hess_struct(fixed_var_func),
                                 num_fixed,
                                 fixed_indices));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_fixed_var_lsq_func_create(SleqpFunc** star,
                                SleqpFunc* func,
                                SleqpParams* params,
                                int num_fixed,
                                const int* fixed_indices,
                                const double* fixed_values)
{
  FixedVarFuncData* func_data;

  assert(sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_LSQ);

  const int num_variables   = sleqp_func_num_vars(func);
  const int num_constraints = sleqp_func_num_cons(func);

  assert(num_fixed >= 0);
  assert(num_fixed <= num_variables);

  SLEQP_CALL(create_fixed_var_func_data(&func_data,
                                        func,
                                        num_fixed,
                                        fixed_indices,
                                        fixed_values));

  SleqpLSQCallbacks callbacks = {.set_value       = fixed_var_func_set,
                                 .lsq_nonzeros    = fixed_lsq_func_nonzeros,
                                 .lsq_residuals   = fixed_lsq_func_residuals,
                                 .lsq_jac_forward = fixed_lsq_func_jac_forward,
                                 .lsq_jac_adjoint = fixed_lsq_func_jac_adjoint,
                                 .cons_val        = fixed_var_cons_val,
                                 .cons_jac        = fixed_var_cons_jac,
                                 .func_free       = fixed_func_free};

  const double levenberg_marquardt
    = sleqp_lsq_func_get_levenberg_marquardt(func);
  const int num_residuals = sleqp_lsq_func_num_residuals(func);

  SLEQP_CALL(sleqp_lsq_func_create(star,
                                   &callbacks,
                                   num_variables - num_fixed,
                                   num_constraints,
                                   num_residuals,
                                   levenberg_marquardt,
                                   params,
                                   (void*)func_data));

  SleqpFunc* fixed_var_func = *star;

  SLEQP_CALL(sleqp_func_flags_add(fixed_var_func, SLEQP_FUNC_INTERNAL));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_fixed_var_dyn_func_create(SleqpFunc** star,
                                SleqpFunc* func,
                                int num_fixed,
                                const int* fixed_indices,
                                const double* fixed_values)
{
  FixedVarFuncData* func_data;

  assert(sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_DYNAMIC);

  const int num_variables   = sleqp_func_num_vars(func);
  const int num_constraints = sleqp_func_num_cons(func);

  assert(num_fixed >= 0);
  assert(num_fixed <= num_variables);

  SLEQP_CALL(create_fixed_var_func_data(&func_data,
                                        func,
                                        num_fixed,
                                        fixed_indices,
                                        fixed_values));

  SleqpDynFuncCallbacks callbacks
    = {.set_value        = fixed_var_func_set,
       .nonzeros         = fixed_var_func_nonzeros,
       .set_error_bound  = fixed_dyn_func_set_error_bound,
       .set_obj_weight   = fixed_dyn_func_set_obj_weight,
       .set_cons_weights = fixed_dyn_func_set_cons_weights,
       .eval             = fixed_dyn_func_eval,
       .obj_grad         = fixed_var_obj_grad,
       .cons_jac         = fixed_var_cons_jac,
       .hess_prod        = fixed_var_hess_prod,
       .func_free        = fixed_func_free};

  SLEQP_CALL(sleqp_dyn_func_create(star,
                                   &callbacks,
                                   num_variables - num_fixed,
                                   num_constraints,
                                   (void*)func_data));

  SleqpFunc* fixed_var_func = *star;

  SLEQP_CALL(sleqp_func_flags_copy(func, fixed_var_func, SLEQP_FUNC_HESS_PSD));

  SLEQP_CALL(sleqp_func_flags_add(fixed_var_func, SLEQP_FUNC_INTERNAL));

  SLEQP_CALL(
    create_fixed_var_hess_struct(sleqp_func_hess_struct(func),
                                 sleqp_func_hess_struct(fixed_var_func),
                                 num_fixed,
                                 fixed_indices));

  return SLEQP_OKAY;
}
