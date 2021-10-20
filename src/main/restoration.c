#include "restoration.h"

#include "cmp.h"
#include "mem.h"
#include "lsq.h"

#include "sparse/sparse_matrix.h"

typedef struct {
  SleqpSparseVec* var_primal;
  SleqpSparseVec* cons_primal;

  SleqpSparseVec* var_forward;
  SleqpSparseVec* cons_forward;

  double* forward_cache;
  SleqpSparseVec* forward_product;

  SleqpSparseVec* adjoint_product;

  SleqpSparseVec* cons_val;
  SleqpSparseMatrix* cons_jac;

  bool has_cons_val;
  bool has_cons_jac;

  SleqpProblem* problem;
  SleqpParams* params;

} FuncData;

static SLEQP_RETCODE
split_primal(FuncData* func_data,
             const SleqpSparseVec* primal,
             SleqpSparseVec* var_primal,
             SleqpSparseVec* cons_primal)
{
  SLEQP_CALL(sleqp_sparse_vector_clear(var_primal));
  SLEQP_CALL(sleqp_sparse_vector_clear(cons_primal));

  SLEQP_CALL(sleqp_sparse_vector_reserve(var_primal,
                                         SLEQP_MIN(primal->nnz, var_primal->dim)));

  SLEQP_CALL(sleqp_sparse_vector_reserve(cons_primal,
                                         SLEQP_MIN(primal->nnz, cons_primal->dim)));

  for(int k = 0; k < primal->nnz; ++k)
  {
    if(primal->indices[k] < var_primal->dim)
    {
      SLEQP_CALL(sleqp_sparse_vector_push(var_primal,
                                          primal->indices[k],
                                          primal->data[k]));
    }
    else
    {
      SLEQP_CALL(sleqp_sparse_vector_push(cons_primal,
                                          primal->indices[k] - var_primal->dim,
                                          primal->data[k]));
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
restoration_func_set(SleqpFunc* func,
                     SleqpSparseVec* value,
                     SLEQP_VALUE_REASON reason,
                     bool* reject,
                     int* func_grad_nnz,
                     int* cons_val_nnz,
                     int* cons_jac_nnz,
                     void* data)
{
  FuncData* func_data = (FuncData*) data;

  SLEQP_CALL(split_primal(func_data,
                          value,
                          func_data->var_primal,
                          func_data->cons_primal));

  func_data->has_cons_val = false;
  func_data->has_cons_jac = false;

  int restoration_func_grad_nnz;
  int restoration_cons_val_nnz;
  int restoration_cons_jac_nnz;

  SLEQP_CALL(sleqp_problem_set_value(func_data->problem,
                                     func_data->var_primal,
                                     reason,
                                     reject,
                                     &restoration_func_grad_nnz,
                                     &restoration_cons_val_nnz,
                                     &restoration_cons_jac_nnz));

  SLEQP_CALL(sleqp_sparse_vector_reserve(func_data->cons_val,
                                         restoration_cons_val_nnz));

  SLEQP_CALL(sleqp_sparse_matrix_reserve(func_data->cons_jac,
                                         restoration_cons_jac_nnz));

  *func_grad_nnz = sleqp_func_get_num_variables(func);
  *cons_val_nnz = 0;
  *cons_jac_nnz = 0;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_cons_val(FuncData* func_data)
{
  if(!func_data->has_cons_val)
  {
    SLEQP_CALL(sleqp_problem_cons_val(func_data->problem,
                                      NULL,
                                      func_data->cons_val));
    func_data->has_cons_val = true;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
restoration_lsq_residuals(SleqpFunc* func,
                          SleqpSparseVec* residual,
                          void* data)
{
  FuncData* func_data =(FuncData*) data;

  SLEQP_CALL(compute_cons_val(func_data));

  const double zero_eps = sleqp_params_get(func_data->params,
                                           SLEQP_PARAM_ZERO_EPS);

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(func_data->cons_val,
                                            func_data->cons_primal,
                                            1.,
                                            -1.,
                                            zero_eps,
                                            residual));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_cons_jac(FuncData* func_data)
{
  if(!func_data->has_cons_jac)
  {
    SLEQP_CALL(sleqp_problem_cons_jac(func_data->problem,
                                      NULL,
                                      func_data->cons_jac));

    func_data->has_cons_jac = true;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
restoration_lsq_jac_forward(SleqpFunc* func,
                            const SleqpSparseVec* forward_direction,
                            SleqpSparseVec* product,
                            void* data)
{
  FuncData* func_data =(FuncData*) data;

  const double zero_eps = sleqp_params_get(func_data->params,
                                           SLEQP_PARAM_ZERO_EPS);

  const int num_constraints = sleqp_problem_num_constraints(func_data->problem);

  SLEQP_CALL(compute_cons_jac(func_data));

  SLEQP_CALL(split_primal(func_data,
                          forward_direction,
                          func_data->var_forward,
                          func_data->cons_forward));

  SLEQP_CALL(sleqp_sparse_matrix_vector_product(func_data->cons_jac,
                                                func_data->var_forward,
                                                func_data->forward_cache));

  SLEQP_CALL(sleqp_sparse_vector_from_raw(func_data->forward_product,
                                          func_data->forward_cache,
                                          num_constraints,
                                          zero_eps));

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(func_data->forward_product,
                                            func_data->cons_forward,
                                            1.,
                                            -1.,
                                            zero_eps,
                                            product));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
concat_adjoint(const SleqpSparseVec* adjoint_product,
               const SleqpSparseVec* adjoint_direction,
               SleqpSparseVec* result)
{
  SLEQP_CALL(sleqp_sparse_vector_clear(result));

  const int total_nnz = adjoint_product->nnz + adjoint_direction->nnz;

  SLEQP_CALL(sleqp_sparse_vector_reserve(result,
                                         SLEQP_MIN(total_nnz, result->dim)));

  for(int k = 0; k < adjoint_product->nnz;++k)
  {
    SLEQP_CALL(sleqp_sparse_vector_push(result,
                                        adjoint_product->indices[k],
                                        adjoint_product->data[k]));
  }

  const int offset = adjoint_product->dim;

  for(int k = 0; k < adjoint_direction->nnz; ++k)
  {
    SLEQP_CALL(sleqp_sparse_vector_push(result,
                                        adjoint_direction->indices[k] + offset,
                                        (-1.) * adjoint_direction->data[k]));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
restoration_lsq_jac_adjoint(SleqpFunc* func,
                            const SleqpSparseVec* adjoint_direction,
                            SleqpSparseVec* product,
                            void* data)
{
  FuncData* func_data =(FuncData*) data;

  SLEQP_CALL(compute_cons_jac(func_data));

  const double zero_eps = sleqp_params_get(func_data->params,
                                           SLEQP_PARAM_ZERO_EPS);

  SLEQP_CALL(sleqp_sparse_matrix_trans_vector_product(func_data->cons_jac,
                                                      adjoint_direction,
                                                      zero_eps,
                                                      func_data->adjoint_product));

  SLEQP_CALL(concat_adjoint(func_data->adjoint_product,
                            adjoint_direction,
                            product));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
restoration_func_free(void* data)
{
  FuncData* func_data = (FuncData*) data;

  SLEQP_CALL(sleqp_params_release(&func_data->params));

  SLEQP_CALL(sleqp_problem_release(&func_data->problem));

  SLEQP_CALL(sleqp_sparse_vector_free(&func_data->adjoint_product));

  SLEQP_CALL(sleqp_sparse_vector_free(&func_data->forward_product));

  sleqp_free(&func_data->forward_cache);

  SLEQP_CALL(sleqp_sparse_vector_free(&func_data->cons_forward));

  SLEQP_CALL(sleqp_sparse_vector_free(&func_data->var_forward));

  SLEQP_CALL(sleqp_sparse_vector_free(&func_data->cons_primal));

  SLEQP_CALL(sleqp_sparse_vector_free(&func_data->var_primal));

  SLEQP_CALL(sleqp_sparse_matrix_release(&func_data->cons_jac));

  SLEQP_CALL(sleqp_sparse_vector_free(&func_data->cons_val));

  sleqp_free(&func_data);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
func_data_create(FuncData** star,
                 SleqpProblem* problem,
                 SleqpParams* params)
{
  SLEQP_CALL(sleqp_malloc(star));

  FuncData* func_data = *star;

  const int num_variables = sleqp_problem_num_variables(problem);
  const int num_constraints = sleqp_problem_num_constraints(problem);

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&func_data->cons_val,
                                              num_constraints));

  SLEQP_CALL(sleqp_sparse_matrix_create(&func_data->cons_jac,
                                        num_constraints,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&func_data->var_primal,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&func_data->cons_primal,
                                              num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&func_data->var_forward,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&func_data->cons_forward,
                                              num_constraints));

  SLEQP_CALL(sleqp_alloc_array(&func_data->forward_cache,
                               num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&func_data->forward_product,
                                              num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&func_data->adjoint_product,
                                              num_variables));

  SLEQP_CALL(sleqp_problem_capture(problem));
  func_data->problem = problem;

  SLEQP_CALL(sleqp_params_capture(params));
  func_data->params = params;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
restoration_func_create(SleqpFunc** star,
                        SleqpParams* params,
                        SleqpProblem* problem)
{
  SleqpLSQCallbacks callbacks = {
    .set_value       = restoration_func_set,
    .lsq_residuals   = restoration_lsq_residuals,
    .lsq_jac_forward = restoration_lsq_jac_forward,
    .lsq_jac_adjoint = restoration_lsq_jac_adjoint,
    .cons_val        = NULL,
    .cons_jac        = NULL,
    .func_free       = restoration_func_free
  };

  const int num_variables = sleqp_problem_num_variables(problem);
  const int num_constraints = sleqp_problem_num_constraints(problem);

  const int restoration_num_variables = num_variables + num_constraints;
  const int restoration_num_constraints = 0;
  const int restoration_num_residuals = num_constraints;

  FuncData* func_data;

  SLEQP_CALL(func_data_create(&func_data, problem, params));

  SLEQP_CALL(sleqp_lsq_func_create(star,
                                   &callbacks,
                                   restoration_num_variables,
                                   restoration_num_constraints,
                                   restoration_num_residuals,
                                   0.,
                                   params,
                                   func_data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_restoration_problem_create(SleqpProblem** star,
                                               SleqpParams* params,
                                               SleqpProblem* problem)
{
  const int num_variables = sleqp_problem_num_variables(problem);
  const int num_constraints = sleqp_problem_num_constraints(problem);

  SleqpFunc* restoration_func;

  SleqpSparseVec* restoration_var_lb;
  SleqpSparseVec* restoration_var_ub;
  SleqpSparseVec* empty;

  const int restoration_num_variables = num_variables + num_constraints;

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&restoration_var_lb,
                                              restoration_num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&restoration_var_ub,
                                              restoration_num_variables));

  SLEQP_CALL(sleqp_sparse_vector_concat(sleqp_problem_var_lb(problem),
                                        sleqp_problem_cons_lb(problem),
                                        restoration_var_lb));

  SLEQP_CALL(sleqp_sparse_vector_concat(sleqp_problem_var_ub(problem),
                                        sleqp_problem_cons_ub(problem),
                                        restoration_var_ub));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&empty, 0));

  SLEQP_CALL(restoration_func_create(&restoration_func, params, problem));

  SLEQP_CALL(sleqp_problem_create_simple(star,
                                         restoration_func,
                                         params,
                                         restoration_var_lb,
                                         restoration_var_ub,
                                         empty,
                                         empty));

  SLEQP_CALL(sleqp_sparse_vector_free(&empty));

  SLEQP_CALL(sleqp_sparse_vector_free(&restoration_var_ub));
  SLEQP_CALL(sleqp_sparse_vector_free(&restoration_var_lb));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_restoration_problem_transform(SleqpProblem* problem,
                                                  const SleqpSparseVec* primal,
                                                  const SleqpSparseVec* cons_val,
                                                  SleqpSparseVec* result)
{
  const int num_variables = sleqp_problem_num_variables(problem);
  const int num_constraints = sleqp_problem_num_constraints(problem);

  assert(primal->dim == num_variables);
  assert(cons_val->dim == num_constraints);

  assert(result->dim == (num_variables + num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_clear(result));

  SLEQP_CALL(sleqp_sparse_vector_reserve(result,
                                         primal->nnz + cons_val->nnz));

  for(int k = 0; k < primal->nnz; ++k)
  {
    SLEQP_CALL(sleqp_sparse_vector_push(result,
                                        primal->indices[k],
                                        primal->data[k]));
  }

  int k_lb = 0, k_c = 0, k_ub = 0;

  const SleqpSparseVec* c = cons_val;
  const SleqpSparseVec* lb = sleqp_problem_cons_lb(problem);
  const SleqpSparseVec* ub = sleqp_problem_cons_ub(problem);

  while(k_c < c->nnz || k_lb < lb->nnz || k_ub < ub->nnz)
  {
    bool valid_c = (k_c < c->nnz);
    bool valid_lb = (k_lb < lb->nnz);
    bool valid_ub = (k_ub < ub->nnz);

    int i_c = valid_c ? c->indices[k_c] : num_constraints + 1;
    int i_lb = valid_lb ? lb->indices[k_lb] : num_constraints + 1;
    int i_ub = valid_ub ? ub->indices[k_ub] : num_constraints + 1;

    int i_combined;

    i_combined = SLEQP_MIN(i_lb, i_ub);
    i_combined = SLEQP_MIN(i_combined, i_c);

    valid_c = valid_c && (i_c == i_combined);
    valid_lb = valid_lb && (i_lb == i_combined);
    valid_ub = valid_ub && (i_ub == i_combined);

    double c_val = valid_c ? c->data[k_c] : 0.;
    double lb_val = valid_lb ? lb->data[k_lb] : 0.;
    double ub_val = valid_ub ? ub->data[k_ub] : 0.;

    c_val = SLEQP_MAX(SLEQP_MIN(c_val, ub_val), lb_val);

    if(c_val != 0.)
    {
      SLEQP_CALL(sleqp_sparse_vector_push(result,
                                          num_variables + i_combined,
                                          c_val));
    }

    if(valid_c)
    {
      ++k_c;
    }

    if(valid_lb)
    {
      ++k_lb;
    }

    if(valid_ub)
    {
      ++k_ub;
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_restoration_problem_restore(SleqpProblem* problem,
                                                const SleqpSparseVec* input,
                                                SleqpSparseVec* result)
{
  const int num_variables = sleqp_problem_num_variables(problem);
  const int num_constraints = sleqp_problem_num_constraints(problem);

  assert(input->dim == (num_variables + num_constraints));
  assert(result->dim == num_variables);

  SLEQP_CALL(sleqp_sparse_vector_clear(result));

  SLEQP_CALL(sleqp_sparse_vector_reserve(result,
                                         SLEQP_MIN(input->dim, num_variables)));

  for(int k = 0; k < input->nnz; ++k)
  {
    const int j = input->indices[k];

    if(j >= num_variables)
    {
      break;
    }

    SLEQP_CALL(sleqp_sparse_vector_push(result,
                                        j,
                                        input->data[k]));
  }

  return SLEQP_OKAY;
}