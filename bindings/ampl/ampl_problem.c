#include "ampl_problem.h"

#include <assert.h>

#include "ampl_mem.h"
#include "ampl_util.h"

#define SLEQP_AMPL_ERROR_CHECK(error)                                          \
  do                                                                           \
  {                                                                            \
    if (error)                                                                 \
    {                                                                          \
      sleqp_log_error("Encountered error during evaluation:");                 \
      return SLEQP_INTERNAL_ERROR;                                             \
    }                                                                          \
  } while (false)

typedef struct AmplFuncData
{
  SleqpAmplData* ampl_data;
  double zero_eps;

  double* x;
  double* cons_vals;
  double* func_grad;

  double* direction;
  double* multipliers;
  double* hessian_product;

  int jac_nnz;

  int* jac_rows;
  int* jac_cols;
  double* jac_vals;

  bool inverted_obj;
  double offset;

} AmplFuncData;

static SLEQP_RETCODE
ampl_func_data_create(AmplFuncData** star,
                      SleqpAmplData* ampl_data,
                      double zero_eps)
{
  ASL* asl = ampl_data->asl;
  SLEQP_CALL(sleqp_ampl_malloc(star));

  AmplFuncData* data = *star;

  *data = (AmplFuncData){0};

  const int num_variables   = ampl_data->num_variables;
  const int num_constraints = ampl_data->num_constraints;

  data->ampl_data = ampl_data;
  data->zero_eps  = zero_eps;
  data->jac_nnz   = nzc;

  data->inverted_obj = sleqp_ampl_max_problem(asl);

  if (n_obj > 0)
  {
    data->offset = objconst(0);
  }
  else
  {
    data->offset = 0.;
  }

  SLEQP_CALL(sleqp_ampl_alloc_array(&data->x, num_variables));
  SLEQP_CALL(sleqp_ampl_alloc_array(&data->cons_vals, num_constraints));
  SLEQP_CALL(sleqp_ampl_alloc_array(&data->multipliers, num_constraints));

  SLEQP_CALL(sleqp_ampl_alloc_array(&data->func_grad, num_variables));

  SLEQP_CALL(sleqp_ampl_alloc_array(&data->direction, num_variables));
  SLEQP_CALL(sleqp_ampl_alloc_array(&data->hessian_product, num_variables));

  SLEQP_CALL(sleqp_ampl_alloc_array(&data->jac_rows, data->jac_nnz));
  SLEQP_CALL(sleqp_ampl_alloc_array(&data->jac_cols, data->jac_nnz));
  SLEQP_CALL(sleqp_ampl_alloc_array(&data->jac_vals, data->jac_nnz));

  for (int i = 0; i < num_constraints; ++i)
  {
    for (cgrad* cg = Cgrad[i]; cg; cg = cg->next)
    {
      data->jac_rows[cg->goff] = i;
      data->jac_cols[cg->goff] = cg->varno;
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
ampl_func_data_free(void* func_data)
{
  AmplFuncData* data  = (AmplFuncData*)func_data;
  AmplFuncData** star = &data;

  sleqp_ampl_free(&data->jac_vals);
  sleqp_ampl_free(&data->jac_cols);
  sleqp_ampl_free(&data->jac_rows);
  sleqp_ampl_free(&data->hessian_product);
  sleqp_ampl_free(&data->direction);

  sleqp_ampl_free(&data->func_grad);
  sleqp_ampl_free(&data->multipliers);
  sleqp_ampl_free(&data->cons_vals);
  sleqp_ampl_free(&data->x);

  sleqp_ampl_free(star);

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
  *cons_jac_nnz = data->jac_nnz;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
ampl_obj_val(SleqpFunc* func, double* func_val, void* func_data)
{
  AmplFuncData* data = (AmplFuncData*)func_data;
  ASL* asl           = data->ampl_data->asl;

  fint nerror = 0;
  *func_val   = objval(0, data->x, &nerror);

  SLEQP_AMPL_ERROR_CHECK(nerror);

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

  fint nerror = 0;
  objgrd(0, data->x, data->func_grad, &nerror);

  SLEQP_AMPL_ERROR_CHECK(nerror);

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
ampl_cons_val(SleqpFunc* func,
              const SleqpSparseVec* cons_indices,
              SleqpSparseVec* cons_val,
              void* func_data)
{
  AmplFuncData* data = (AmplFuncData*)func_data;
  ASL* asl           = data->ampl_data->asl;

  fint nerror = 0;
  conval(data->x, data->cons_vals, &nerror);

  SLEQP_AMPL_ERROR_CHECK(nerror);

  SLEQP_CALL(sleqp_sparse_vector_from_raw(cons_val,
                                          data->cons_vals,
                                          data->ampl_data->num_constraints,
                                          data->zero_eps));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
ampl_cons_jac(SleqpFunc* func,
              const SleqpSparseVec* cons_indices,
              SleqpSparseMatrix* cons_jac,
              void* func_data)
{
  AmplFuncData* data = (AmplFuncData*)func_data;
  ASL* asl           = data->ampl_data->asl;

  fint nerror = 0;
  jacval(data->x, data->jac_vals, &nerror);

  SLEQP_AMPL_ERROR_CHECK(nerror);

  int last_col = 0;

  SLEQP_CALL(sleqp_sparse_matrix_reserve(cons_jac, data->jac_nnz));

  for (int i = 0; i < data->jac_nnz; ++i)
  {
    int row    = data->jac_rows[i];
    int col    = data->jac_cols[i];
    double val = data->jac_vals[i];

    while (col > last_col)
    {
      SLEQP_CALL(sleqp_sparse_matrix_push_column(cons_jac, ++last_col));
    }

    last_col = col;

    SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac, row, col, val));
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
                       SleqpParams* params)
{
  AmplFuncData* data;

  const int num_variables   = ampl_data->num_variables;
  const int num_constraints = ampl_data->num_constraints;

  const double zero_eps = sleqp_params_value(params, SLEQP_PARAM_ZERO_EPS);

  SLEQP_CALL(ampl_func_data_create(&data, ampl_data, zero_eps));

  SleqpFuncCallbacks callbacks
    = {.set_value = ampl_func_set,
       .obj_val   = ampl_obj_val,
       .obj_grad  = ampl_obj_grad,
       .cons_val  = ampl_data->is_constrained ? ampl_cons_val : NULL,
       .cons_jac  = ampl_data->is_constrained ? ampl_cons_jac : NULL,
       .hess_prod = ampl_func_hess_product,
       .func_free = ampl_func_data_free};

  SLEQP_CALL(
    sleqp_func_create(star, &callbacks, num_variables, num_constraints, data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_ampl_problem_create(SleqpProblem** star,
                          SleqpAmplData* data,
                          FILE* nl,
                          SleqpParams* params)
{
  const int num_variables   = data->num_variables;
  const int num_constraints = data->num_constraints;

  ASL* asl = data->asl;

  // set ASL pointer to allocated data
  X0    = data->x;
  LUv   = data->var_lb;
  Uvx   = data->var_ub;
  LUrhs = data->cons_lb;
  Urhsx = data->cons_ub;

  // read remaining data from stub and set functions
  int retcode = pfgh_read(nl, ASL_return_read_err);
  if (retcode != ASL_readerr_none)
  {
    sleqp_log_error("Error in reading nl file %d", retcode);
  }

  SLEQP_CALL(map_ampl_inf(data->var_lb, num_variables));
  SLEQP_CALL(map_ampl_inf(data->var_ub, num_variables));

  SLEQP_CALL(map_ampl_inf(data->cons_lb, num_constraints));
  SLEQP_CALL(map_ampl_inf(data->cons_ub, num_constraints));

  const double zero_eps = sleqp_params_value(params, SLEQP_PARAM_ZERO_EPS);

  SleqpSparseVec* var_lb;
  SleqpSparseVec* var_ub;

  SleqpSparseVec* cons_lb;
  SleqpSparseVec* cons_ub;

  SleqpFunc* func;

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&var_lb, num_variables));
  SLEQP_CALL(sleqp_sparse_vector_create_empty(&var_ub, num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&cons_lb, num_constraints));
  SLEQP_CALL(sleqp_sparse_vector_create_empty(&cons_ub, num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_from_raw(var_lb,
                                          data->var_lb,
                                          num_variables,
                                          zero_eps));
  SLEQP_CALL(sleqp_sparse_vector_from_raw(var_ub,
                                          data->var_ub,
                                          num_variables,
                                          zero_eps));

  SLEQP_CALL(sleqp_sparse_vector_from_raw(cons_lb,
                                          data->cons_lb,
                                          num_constraints,
                                          zero_eps));
  SLEQP_CALL(sleqp_sparse_vector_from_raw(cons_ub,
                                          data->cons_ub,
                                          num_constraints,
                                          zero_eps));

  SLEQP_CALL(sleqp_ampl_func_create(&func, data, params));

  SLEQP_CALL(sleqp_problem_create_simple(star,
                                         func,
                                         params,
                                         var_lb,
                                         var_ub,
                                         cons_lb,
                                         cons_ub));

  SLEQP_CALL(sleqp_func_release(&func));

  SLEQP_CALL(sleqp_sparse_vector_free(&cons_ub));
  SLEQP_CALL(sleqp_sparse_vector_free(&cons_lb));

  SLEQP_CALL(sleqp_sparse_vector_free(&var_ub));
  SLEQP_CALL(sleqp_sparse_vector_free(&var_lb));

  return SLEQP_OKAY;
}
