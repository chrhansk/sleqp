#include "sleqp_fixed_var_func.h"

#include "sleqp_assert.h"
#include "sleqp_cmp.h"
#include "sleqp_mem.h"

#include "preprocessor/sleqp_preprocessing.h"

typedef struct FixedVarFuncData
{
  int num_fixed;
  double* fixed_values;
  int* fixed_indices;

  SleqpFunc* func;

  SleqpSparseVec* values;

  SleqpSparseVec* grad;

  SleqpSparseVec* direction;
  SleqpSparseVec* product;

  SleqpSparseMatrix* jacobian;

} FixedVarFuncData;


static
SLEQP_RETCODE fixed_var_func_set(SleqpFunc* func,
                                 SleqpSparseVec* value,
                                 SLEQP_VALUE_REASON reason,
                                 int* func_grad_nnz,
                                 int* cons_val_nnz,
                                 int* cons_jac_nnz,
                                 void* data)
{
  FixedVarFuncData* func_data = (FixedVarFuncData*) data;

  SLEQP_CALL(sleqp_preprocessing_merge_fixed_entries(value, func_data->values,
                                                     func_data->num_fixed,
                                                     func_data->fixed_indices,
                                                     func_data->fixed_values));

  SLEQP_CALL(sleqp_func_set_value(func_data->func,
                                  func_data->values,
                                  reason,
                                  func_grad_nnz,
                                  cons_val_nnz,
                                  cons_jac_nnz));

  SLEQP_CALL(sleqp_sparse_vector_reserve(func_data->grad,
                                         *func_grad_nnz));

  SLEQP_CALL(sleqp_sparse_matrix_reserve(func_data->jacobian,
                                         *cons_jac_nnz));

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE fixed_var_func_val(SleqpFunc* func,
                                 double* func_val,
                                 void* data)
{
  FixedVarFuncData* func_data = (FixedVarFuncData*) data;

  return sleqp_func_val(func_data->func, func_val);
}

static
SLEQP_RETCODE fixed_var_func_grad(SleqpFunc* func,
                                  SleqpSparseVec* func_grad,
                                  void* data)
{
  FixedVarFuncData* func_data = (FixedVarFuncData*) data;

  SLEQP_CALL(sleqp_func_grad(func_data->func, func_data->grad));

  SLEQP_CALL(sleqp_preprocessing_remove_fixed_entries(func_data->grad,
                                                      func_grad,
                                                      func_data->num_fixed,
                                                      func_data->fixed_indices));

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE fixed_var_cons_val(SleqpFunc* func,
                                 const SleqpSparseVec* cons_indices,
                                 SleqpSparseVec* cons_val,
                                 void* data)
{
  FixedVarFuncData* func_data = (FixedVarFuncData*) data;

  return sleqp_func_cons_val(func_data->func, cons_indices, cons_val);
}

static
SLEQP_RETCODE fixed_var_cons_jac(SleqpFunc* func,
                                 const SleqpSparseVec* cons_indices,
                                 SleqpSparseMatrix* cons_jac,
                                 void* data)
{
  FixedVarFuncData* func_data = (FixedVarFuncData*) data;

  SLEQP_CALL(sleqp_func_cons_jac(func_data->func,
                                 cons_indices,
                                 func_data->jacobian));

  SLEQP_CALL(sleqp_preprocessing_remove_fixed_matrix_entries(func_data->jacobian,
                                                             cons_jac,
                                                             func_data->num_fixed,
                                                             func_data->fixed_indices));

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE fixed_var_hess_prod(SleqpFunc* func,
                                  const double* func_dual,
                                  const SleqpSparseVec* direction,
                                  const SleqpSparseVec* cons_duals,
                                  SleqpSparseVec* product,
                                  void* data)
{
  FixedVarFuncData* func_data = (FixedVarFuncData*) data;

  SLEQP_CALL(sleqp_preprocessing_add_zero_entries(direction,
                                                  func_data->direction,
                                                  func_data->num_fixed,
                                                  func_data->fixed_indices));

  SLEQP_CALL(sleqp_func_hess_prod(func_data->func,
                                  func_dual,
                                  func_data->direction,
                                  cons_duals,
                                  func_data->product));

  SLEQP_CALL(sleqp_preprocessing_remove_fixed_entries(func_data->product,
                                                      product,
                                                      func_data->num_fixed,
                                                      func_data->fixed_indices));

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE fixed_func_free(void* data)
{
  FixedVarFuncData* func_data = (FixedVarFuncData*) data;

  SLEQP_CALL(sleqp_sparse_matrix_release(&func_data->jacobian));

  SLEQP_CALL(sleqp_sparse_vector_free(&func_data->product));

  SLEQP_CALL(sleqp_sparse_vector_free(&func_data->direction));

  SLEQP_CALL(sleqp_sparse_vector_free(&func_data->grad));

  SLEQP_CALL(sleqp_sparse_vector_free(&func_data->values));

  SLEQP_CALL(sleqp_func_release(&func_data->func));

  sleqp_free(&func_data->fixed_values);

  sleqp_free(&func_data->fixed_indices);

  sleqp_free(&func_data);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_fixed_var_func_create(SleqpFunc** star,
                                          SleqpFunc* func,
                                          int num_fixed,
                                          const int* fixed_indices,
                                          const double* fixed_values)
{
  FixedVarFuncData* func_data;

  const int num_variables = sleqp_func_get_num_variables(func);
  const int num_constraints = sleqp_func_get_num_constraints(func);

  assert(num_fixed >= 0);
  assert(num_fixed <= num_variables);

  SLEQP_CALL(sleqp_malloc(&func_data));

  *func_data = (FixedVarFuncData) {0};

  func_data->num_fixed = num_fixed;

  SLEQP_CALL(sleqp_alloc_array(&func_data->fixed_indices,
                               num_fixed));

  SLEQP_CALL(sleqp_alloc_array(&func_data->fixed_values,
                               num_fixed));

  int last_index = -1;

  for(int j = 0; j < num_fixed; ++j)
  {
    const int fixed_index = fixed_indices[j];

    sleqp_assert(fixed_index >= 0);
    sleqp_assert(fixed_index < num_variables);
    sleqp_assert(fixed_index > last_index);

    func_data->fixed_indices[j] = fixed_index;

    last_index = fixed_index;
  }

  for(int j = 0; j < num_fixed; ++j)
  {
    const double fixed_value = fixed_values[j];

    sleqp_assert(sleqp_is_finite(fixed_value));

    func_data->fixed_values[j] = fixed_value;
  }

  func_data->func = func;

  SLEQP_CALL(sleqp_func_capture(func));

  SLEQP_CALL(sleqp_sparse_vector_create_full(&func_data->values,
                                             num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_full(&func_data->grad,
                                             num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_full(&func_data->direction,
                                             num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_full(&func_data->product,
                                             num_variables));

  SLEQP_CALL(sleqp_sparse_matrix_create(&func_data->jacobian,
                                        num_constraints,
                                        num_variables,
                                        0));

  SleqpFuncCallbacks callbacks = {
    .set_value = fixed_var_func_set,
    .func_val = fixed_var_func_val,
    .func_grad = fixed_var_func_grad,
    .cons_val = fixed_var_cons_val,
    .cons_jac = fixed_var_cons_jac,
    .hess_prod = fixed_var_hess_prod,
    .func_free = fixed_func_free
  };

  SLEQP_CALL(sleqp_func_create(star,
                               &callbacks,
                               num_variables - num_fixed,
                               num_constraints,
                               (void*) func_data));

  return SLEQP_OKAY;
}