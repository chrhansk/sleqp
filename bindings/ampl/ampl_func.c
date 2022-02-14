#include "ampl_func.h"

#include "ampl_util.h"

#define SLEQP_AMPL_ERROR_CHECK(errorptr)                                       \
  do                                                                           \
  {                                                                            \
    if (errorptr && (*errorptr != 0))                                          \
    {                                                                          \
      sleqp_raise(SLEQP_INTERNAL_ERROR,                                        \
                  "Error during evaluation. "                                  \
                  "Run with \"halt_on_ampl_error yes\" to see details.");      \
    }                                                                          \
  } while (false)

typedef struct AmplFuncData
{
  SleqpAmplData* ampl_data;
  double zero_eps;

  double* x;
  double* func_grad;

  double* direction;
  double* multipliers;
  double* hessian_product;

  bool inverted_obj;
  double offset;

  fint error;
  fint* nerror;

} AmplFuncData;

static SLEQP_RETCODE
ampl_func_data_create(AmplFuncData** star,
                      SleqpAmplData* ampl_data,
                      double zero_eps,
                      bool halt_on_error)
{
  ASL* asl = ampl_data->asl;
  SLEQP_CALL(sleqp_malloc(star));

  AmplFuncData* data = *star;

  *data = (AmplFuncData){0};

  const int num_variables = ampl_data->num_variables;
  const int num_general   = ampl_data->num_constraints - ampl_data->num_linear;

  data->ampl_data = ampl_data;
  data->zero_eps  = zero_eps;

  data->inverted_obj = sleqp_ampl_max_problem(asl);

  if (halt_on_error)
  {
    data->nerror = NULL;
  }
  else
  {
    data->nerror = &data->error;
  }

  if (n_obj > 0)
  {
    data->offset = objconst(0);
  }
  else
  {
    data->offset = 0.;
  }

  SLEQP_CALL(sleqp_alloc_array(&data->x, num_variables));
  SLEQP_CALL(sleqp_alloc_array(&data->multipliers, num_general));

  SLEQP_CALL(sleqp_alloc_array(&data->func_grad, num_variables));

  SLEQP_CALL(sleqp_alloc_array(&data->direction, num_variables));
  SLEQP_CALL(sleqp_alloc_array(&data->hessian_product, num_variables));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
ampl_func_data_free(void* func_data)
{
  AmplFuncData* data  = (AmplFuncData*)func_data;
  AmplFuncData** star = &data;

  sleqp_free(&data->hessian_product);
  sleqp_free(&data->direction);

  sleqp_free(&data->func_grad);
  sleqp_free(&data->multipliers);
  sleqp_free(&data->x);

  sleqp_free(star);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
ampl_func_set(SleqpFunc* func,
              SleqpSparseVec* x,
              SLEQP_VALUE_REASON reason,
              bool* reject,
              int* func_grad_nnz,
              int* cons_val_nnz,
              int* cons_jac_nnz,
              void* func_data)
{
  AmplFuncData* data = (AmplFuncData*)func_data;

  SLEQP_CALL(sleqp_sparse_vector_to_raw(x, data->x));

  *func_grad_nnz = data->ampl_data->num_variables;

  *cons_val_nnz = data->ampl_data->num_constraints;
  *cons_jac_nnz = data->ampl_data->jac_nnz;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
ampl_obj_val(SleqpFunc* func, double* func_val, void* func_data)
{
  AmplFuncData* data = (AmplFuncData*)func_data;
  ASL* asl           = data->ampl_data->asl;

  *func_val = objval(0, data->x, data->nerror);

  SLEQP_AMPL_ERROR_CHECK(data->nerror);

  *func_val += data->offset;

  if (data->inverted_obj)
  {
    *func_val *= -1.;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
ampl_obj_grad(SleqpFunc* func, SleqpSparseVec* func_grad, void* func_data)
{
  AmplFuncData* data = (AmplFuncData*)func_data;
  ASL* asl           = data->ampl_data->asl;

  objgrd(0, data->x, data->func_grad, data->nerror);

  SLEQP_AMPL_ERROR_CHECK(data->nerror);

  SLEQP_CALL(sleqp_sparse_vector_from_raw(func_grad,
                                          data->func_grad,
                                          data->ampl_data->num_variables,
                                          data->zero_eps));

  if (data->inverted_obj)
  {
    SLEQP_CALL(sleqp_sparse_vector_scale(func_grad, -1.));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
ampl_cons_val(SleqpFunc* func, SleqpSparseVec* cons_val, void* func_data)
{
  AmplFuncData* data       = (AmplFuncData*)func_data;
  SleqpAmplData* ampl_data = data->ampl_data;
  ASL* asl                 = data->ampl_data->asl;

  conval(data->x, ampl_data->cons_val, data->nerror);

  SLEQP_AMPL_ERROR_CHECK(data->nerror);

  const int num_general = ampl_data->num_constraints - ampl_data->num_linear;

  SLEQP_CALL(sleqp_sparse_vector_from_raw(cons_val,
                                          ampl_data->cons_val,
                                          num_general,
                                          data->zero_eps));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
ampl_cons_jac(SleqpFunc* func, SleqpSparseMatrix* cons_jac, void* func_data)
{
  AmplFuncData* data       = (AmplFuncData*)func_data;
  ASL* asl                 = data->ampl_data->asl;
  SleqpAmplData* ampl_data = data->ampl_data;

  jacval(data->x, ampl_data->jac_vals, data->nerror);

  SLEQP_AMPL_ERROR_CHECK(data->nerror);

  int next_col = 0;

  SLEQP_CALL(sleqp_sparse_matrix_reserve(cons_jac, ampl_data->jac_nnz));

  const int num_general = ampl_data->num_constraints - ampl_data->num_linear;

  for (int i = 0; i < ampl_data->jac_nnz; ++i)
  {
    const int row    = ampl_data->jac_rows[i];
    const int col    = ampl_data->jac_cols[i];
    const double val = ampl_data->jac_vals[i];

    while (col >= next_col)
    {
      SLEQP_CALL(sleqp_sparse_matrix_push_column(cons_jac, next_col++));
    }

    if (row < num_general)
    {
      SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac, row, col, val));
    }
  }

  const int num_cols = data->ampl_data->num_variables;

  while (num_cols > next_col)
  {
    SLEQP_CALL(sleqp_sparse_matrix_push_column(cons_jac, next_col++));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
ampl_func_hess_product(SleqpFunc* func,
                       const double* obj_dual,
                       const SleqpSparseVec* direction,
                       const SleqpSparseVec* cons_duals,
                       SleqpSparseVec* product,
                       void* func_data)
{
  AmplFuncData* data = (AmplFuncData*)func_data;
  ASL* asl           = data->ampl_data->asl;

  SLEQP_CALL(sleqp_sparse_vector_to_raw(direction, data->direction));
  SLEQP_CALL(sleqp_sparse_vector_to_raw(cons_duals, data->multipliers));

  if (data->inverted_obj)
  {
    const int num_cons = sleqp_func_num_cons(func);

    for (int i = 0; i < num_cons; ++i)
    {
      data->multipliers[i] *= -1.;
    }
  }

  double objective_dual = obj_dual ? (*obj_dual) : 0.;

  hvcomp(data->hessian_product,
         data->direction,
         0,
         &objective_dual,
         data->multipliers);

  SLEQP_CALL(sleqp_sparse_vector_from_raw(product,
                                          data->hessian_product,
                                          data->ampl_data->num_variables,
                                          data->zero_eps));

  if (data->inverted_obj)
  {
    SLEQP_CALL(sleqp_sparse_vector_scale(product, -1.));
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_ampl_func_create(SleqpFunc** star,
                       SleqpAmplData* ampl_data,
                       SleqpParams* params,
                       bool halt_on_error)
{
  AmplFuncData* data;

  const int num_variables = ampl_data->num_variables;
  const int num_general   = ampl_data->num_constraints - ampl_data->num_linear;

  const double zero_eps = sleqp_params_value(params, SLEQP_PARAM_ZERO_EPS);

  SLEQP_CALL(ampl_func_data_create(&data, ampl_data, zero_eps, halt_on_error));

  SleqpFuncCallbacks callbacks
    = {.set_value = ampl_func_set,
       .obj_val   = ampl_obj_val,
       .obj_grad  = ampl_obj_grad,
       .cons_val  = ampl_data->is_constrained ? ampl_cons_val : NULL,
       .cons_jac  = ampl_data->is_constrained ? ampl_cons_jac : NULL,
       .hess_prod = ampl_func_hess_product,
       .func_free = ampl_func_data_free};

  SLEQP_CALL(
    sleqp_func_create(star, &callbacks, num_variables, num_general, data));

  return SLEQP_OKAY;
}
