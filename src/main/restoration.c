#include "restoration.h"

#include "cmp.h"
#include "lsq.h"
#include "mem.h"

#include "pub_settings.h"
#include "sparse/mat.h"

typedef struct
{
  SleqpVec* var_primal;
  SleqpVec* cons_primal;

  SleqpVec* var_forward;
  SleqpVec* cons_forward;

  double* forward_cache;
  SleqpVec* forward_product;

  SleqpVec* adjoint_product;

  SleqpVec* cons_val;
  SleqpMat* cons_jac;

  bool has_cons_val;
  bool has_cons_jac;

  SleqpProblem* problem;
  SleqpSettings* settings;

} FuncData;

static SLEQP_RETCODE
split_primal(FuncData* func_data,
             const SleqpVec* primal,
             SleqpVec* var_primal,
             SleqpVec* cons_primal)
{
  SLEQP_CALL(sleqp_vec_clear(var_primal));
  SLEQP_CALL(sleqp_vec_clear(cons_primal));

  SLEQP_CALL(
    sleqp_vec_reserve(var_primal, SLEQP_MIN(primal->nnz, var_primal->dim)));

  SLEQP_CALL(
    sleqp_vec_reserve(cons_primal, SLEQP_MIN(primal->nnz, cons_primal->dim)));

  for (int k = 0; k < primal->nnz; ++k)
  {
    if (primal->indices[k] < var_primal->dim)
    {
      SLEQP_CALL(
        sleqp_vec_push(var_primal, primal->indices[k], primal->data[k]));
    }
    else
    {
      SLEQP_CALL(sleqp_vec_push(cons_primal,
                                primal->indices[k] - var_primal->dim,
                                primal->data[k]));
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
restoration_func_set(SleqpFunc* func,
                     SleqpVec* value,
                     SLEQP_VALUE_REASON reason,
                     bool* reject,
                     void* data)
{
  FuncData* func_data = (FuncData*)data;

  SLEQP_CALL(split_primal(func_data,
                          value,
                          func_data->var_primal,
                          func_data->cons_primal));

  func_data->has_cons_val = false;
  func_data->has_cons_jac = false;

  SLEQP_CALL(sleqp_problem_set_value(func_data->problem,
                                     func_data->var_primal,
                                     reason,
                                     reject));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
restoration_func_nonzeros(SleqpFunc* func,
                          int* residual_nnz,
                          int* jac_fwd_nnz,
                          int* jac_adj_nnz,
                          int* cons_val_nnz,
                          int* cons_jac_nnz,
                          void* data)
{
  *residual_nnz = SLEQP_NONE;
  *jac_fwd_nnz  = SLEQP_NONE;
  *jac_adj_nnz  = SLEQP_NONE;

  FuncData* func_data = (FuncData*)data;

  int problem_obj_grad_nnz  = SLEQP_NONE;
  int problem_cons_val_nnz  = SLEQP_NONE;
  int problem_cons_jac_nnz  = SLEQP_NONE;
  int problem_hess_prod_nnz = SLEQP_NONE;

  SLEQP_CALL(sleqp_problem_nonzeros(func_data->problem,
                                    &problem_obj_grad_nnz,
                                    &problem_cons_val_nnz,
                                    &problem_cons_jac_nnz,
                                    &problem_hess_prod_nnz));

  if (problem_cons_val_nnz != SLEQP_NONE)
  {
    *residual_nnz = problem_cons_val_nnz;

    SLEQP_CALL(sleqp_vec_reserve(func_data->cons_val, problem_cons_val_nnz));
  }

  if (problem_cons_jac_nnz != SLEQP_NONE)
  {
    SLEQP_CALL(sleqp_mat_reserve(func_data->cons_jac, problem_cons_jac_nnz));
  }

  *cons_val_nnz = 0;
  *cons_jac_nnz = 0;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_cons_val(FuncData* func_data)
{
  if (!func_data->has_cons_val)
  {
    SLEQP_CALL(sleqp_problem_cons_val(func_data->problem, func_data->cons_val));
    func_data->has_cons_val = true;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
restoration_lsq_residuals(SleqpFunc* func, SleqpVec* residual, void* data)
{
  FuncData* func_data = (FuncData*)data;

  SLEQP_CALL(compute_cons_val(func_data));

  const double zero_eps
    = sleqp_settings_real_value(func_data->settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  SLEQP_CALL(sleqp_vec_add_scaled(func_data->cons_val,
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
  if (!func_data->has_cons_jac)
  {
    SLEQP_CALL(sleqp_problem_cons_jac(func_data->problem, func_data->cons_jac));

    func_data->has_cons_jac = true;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
restoration_lsq_jac_forward(SleqpFunc* func,
                            const SleqpVec* forward_direction,
                            SleqpVec* product,
                            void* data)
{
  FuncData* func_data = (FuncData*)data;

  const double zero_eps
    = sleqp_settings_real_value(func_data->settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  const int num_constraints = sleqp_problem_num_cons(func_data->problem);

  SLEQP_CALL(compute_cons_jac(func_data));

  SLEQP_CALL(split_primal(func_data,
                          forward_direction,
                          func_data->var_forward,
                          func_data->cons_forward));

  SLEQP_CALL(sleqp_mat_mult_vec(func_data->cons_jac,
                                func_data->var_forward,
                                func_data->forward_cache));

  SLEQP_CALL(sleqp_vec_set_from_raw(func_data->forward_product,
                                    func_data->forward_cache,
                                    num_constraints,
                                    zero_eps));

  SLEQP_CALL(sleqp_vec_add_scaled(func_data->forward_product,
                                  func_data->cons_forward,
                                  1.,
                                  -1.,
                                  zero_eps,
                                  product));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
concat_adjoint(const SleqpVec* adjoint_product,
               const SleqpVec* adjoint_direction,
               SleqpVec* result)
{
  SLEQP_CALL(sleqp_vec_clear(result));

  const int total_nnz = adjoint_product->nnz + adjoint_direction->nnz;

  SLEQP_CALL(sleqp_vec_reserve(result, SLEQP_MIN(total_nnz, result->dim)));

  for (int k = 0; k < adjoint_product->nnz; ++k)
  {
    SLEQP_CALL(sleqp_vec_push(result,
                              adjoint_product->indices[k],
                              adjoint_product->data[k]));
  }

  const int offset = adjoint_product->dim;

  for (int k = 0; k < adjoint_direction->nnz; ++k)
  {
    SLEQP_CALL(sleqp_vec_push(result,
                              adjoint_direction->indices[k] + offset,
                              (-1.) * adjoint_direction->data[k]));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
restoration_lsq_jac_adjoint(SleqpFunc* func,
                            const SleqpVec* adjoint_direction,
                            SleqpVec* product,
                            void* data)
{
  FuncData* func_data = (FuncData*)data;

  SLEQP_CALL(compute_cons_jac(func_data));

  const double zero_eps
    = sleqp_settings_real_value(func_data->settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  SLEQP_CALL(sleqp_mat_mult_vec_trans(func_data->cons_jac,
                                      adjoint_direction,
                                      zero_eps,
                                      func_data->adjoint_product));

  SLEQP_CALL(
    concat_adjoint(func_data->adjoint_product, adjoint_direction, product));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
restoration_func_free(void* data)
{
  FuncData* func_data = (FuncData*)data;

  SLEQP_CALL(sleqp_settings_release(&func_data->settings));

  SLEQP_CALL(sleqp_problem_release(&func_data->problem));

  SLEQP_CALL(sleqp_vec_free(&func_data->adjoint_product));

  SLEQP_CALL(sleqp_vec_free(&func_data->forward_product));

  sleqp_free(&func_data->forward_cache);

  SLEQP_CALL(sleqp_vec_free(&func_data->cons_forward));

  SLEQP_CALL(sleqp_vec_free(&func_data->var_forward));

  SLEQP_CALL(sleqp_vec_free(&func_data->cons_primal));

  SLEQP_CALL(sleqp_vec_free(&func_data->var_primal));

  SLEQP_CALL(sleqp_mat_release(&func_data->cons_jac));

  SLEQP_CALL(sleqp_vec_free(&func_data->cons_val));

  sleqp_free(&func_data);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
func_data_create(FuncData** star, SleqpProblem* problem, SleqpSettings* settings)
{
  SLEQP_CALL(sleqp_malloc(star));

  FuncData* func_data = *star;

  const int num_variables   = sleqp_problem_num_vars(problem);
  const int num_constraints = sleqp_problem_num_cons(problem);

  SLEQP_CALL(sleqp_vec_create_empty(&func_data->cons_val, num_constraints));

  SLEQP_CALL(
    sleqp_mat_create(&func_data->cons_jac, num_constraints, num_variables, 0));

  SLEQP_CALL(sleqp_vec_create_empty(&func_data->var_primal, num_variables));

  SLEQP_CALL(sleqp_vec_create_empty(&func_data->cons_primal, num_constraints));

  SLEQP_CALL(sleqp_vec_create_empty(&func_data->var_forward, num_variables));

  SLEQP_CALL(sleqp_vec_create_empty(&func_data->cons_forward, num_constraints));

  SLEQP_CALL(sleqp_alloc_array(&func_data->forward_cache, num_constraints));

  SLEQP_CALL(
    sleqp_vec_create_empty(&func_data->forward_product, num_constraints));

  SLEQP_CALL(
    sleqp_vec_create_empty(&func_data->adjoint_product, num_variables));

  SLEQP_CALL(sleqp_problem_capture(problem));
  func_data->problem = problem;

  SLEQP_CALL(sleqp_settings_capture(settings));
  func_data->settings = settings;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
restoration_func_create(SleqpFunc** star,
                        SleqpSettings* settings,
                        SleqpProblem* problem)
{
  SleqpLSQCallbacks callbacks = {.set_value       = restoration_func_set,
                                 .lsq_nonzeros    = restoration_func_nonzeros,
                                 .lsq_residuals   = restoration_lsq_residuals,
                                 .lsq_jac_forward = restoration_lsq_jac_forward,
                                 .lsq_jac_adjoint = restoration_lsq_jac_adjoint,
                                 .cons_val        = NULL,
                                 .cons_jac        = NULL,
                                 .func_free       = restoration_func_free};

  const int num_variables   = sleqp_problem_num_vars(problem);
  const int num_constraints = sleqp_problem_num_cons(problem);

  const int restoration_num_variables   = num_variables + num_constraints;
  const int restoration_num_constraints = 0;
  const int restoration_num_residuals   = num_constraints;

  FuncData* func_data;

  SLEQP_CALL(func_data_create(&func_data, problem, settings));

  SLEQP_CALL(sleqp_lsq_func_create(star,
                                   &callbacks,
                                   restoration_num_variables,
                                   restoration_num_constraints,
                                   restoration_num_residuals,
                                   0.,
                                   settings,
                                   func_data));

  SleqpFunc* func = *star;

  SLEQP_CALL(sleqp_func_flags_add(func, SLEQP_FUNC_INTERNAL));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_restoration_func_cons_val(SleqpFunc* restoration_func, SleqpVec** star)
{
  FuncData* func_data = (FuncData*)sleqp_lsq_func_get_data(restoration_func);

  SLEQP_CALL(compute_cons_val(func_data));

  (*star) = func_data->cons_val;

  return SLEQP_OKAY;
}

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_restoration_func_cons_jac(SleqpFunc* restoration_func, SleqpMat** star)
{
  FuncData* func_data = (FuncData*)sleqp_lsq_func_get_data(restoration_func);

  SLEQP_CALL(compute_cons_jac(func_data));

  (*star) = func_data->cons_jac;

  return SLEQP_OKAY;
}

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_restoration_func_init(SleqpFunc* restoration_func,
                            SleqpVec* restoration_primal,
                            SleqpVec* orig_cons_val,
                            SleqpMat* orig_cons_jac)
{
  FuncData* func_data = (FuncData*)sleqp_lsq_func_get_data(restoration_func);

  SLEQP_CALL(split_primal(func_data,
                          restoration_primal,
                          func_data->var_primal,
                          func_data->cons_primal));

  {
    SLEQP_CALL(sleqp_vec_copy(orig_cons_val, func_data->cons_val));
    func_data->has_cons_val = true;
  }

  {
    SLEQP_CALL(sleqp_mat_copy(orig_cons_jac, func_data->cons_jac));
    func_data->has_cons_jac = true;
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_restoration_problem_create(SleqpProblem** star,
                                 SleqpSettings* settings,
                                 SleqpProblem* problem)
{
  const int num_variables   = sleqp_problem_num_vars(problem);
  const int num_constraints = sleqp_problem_num_cons(problem);

  SleqpFunc* restoration_func;

  SleqpVec* restoration_var_lb;
  SleqpVec* restoration_var_ub;
  SleqpVec* empty;

  const int restoration_num_variables = num_variables + num_constraints;

  SLEQP_CALL(
    sleqp_vec_create_empty(&restoration_var_lb, restoration_num_variables));

  SLEQP_CALL(
    sleqp_vec_create_empty(&restoration_var_ub, restoration_num_variables));

  SLEQP_CALL(sleqp_vec_concat(sleqp_problem_vars_lb(problem),
                              sleqp_problem_cons_lb(problem),
                              restoration_var_lb));

  SLEQP_CALL(sleqp_vec_concat(sleqp_problem_vars_ub(problem),
                              sleqp_problem_cons_ub(problem),
                              restoration_var_ub));

  SLEQP_CALL(sleqp_vec_create_empty(&empty, 0));

  SLEQP_CALL(restoration_func_create(&restoration_func, settings, problem));

  SLEQP_CALL(sleqp_problem_create_simple(star,
                                         restoration_func,
                                         settings,
                                         restoration_var_lb,
                                         restoration_var_ub,
                                         empty,
                                         empty));

  SLEQP_CALL(sleqp_func_release(&restoration_func));

  SLEQP_CALL(sleqp_vec_free(&empty));

  SLEQP_CALL(sleqp_vec_free(&restoration_var_ub));
  SLEQP_CALL(sleqp_vec_free(&restoration_var_lb));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_restoration_problem_transform(SleqpProblem* problem,
                                    const SleqpVec* primal,
                                    const SleqpVec* cons_val,
                                    SleqpVec* result)
{
  const int num_variables   = sleqp_problem_num_vars(problem);
  const int num_constraints = sleqp_problem_num_cons(problem);

  assert(primal->dim == num_variables);
  assert(cons_val->dim == num_constraints);

  assert(result->dim == (num_variables + num_constraints));

  SLEQP_CALL(sleqp_vec_clear(result));

  SLEQP_CALL(sleqp_vec_reserve(result, primal->nnz + cons_val->nnz));

  for (int k = 0; k < primal->nnz; ++k)
  {
    SLEQP_CALL(sleqp_vec_push(result, primal->indices[k], primal->data[k]));
  }

  int k_lb = 0, k_c = 0, k_ub = 0;

  const SleqpVec* c  = cons_val;
  const SleqpVec* lb = sleqp_problem_cons_lb(problem);
  const SleqpVec* ub = sleqp_problem_cons_ub(problem);

  while (k_c < c->nnz || k_lb < lb->nnz || k_ub < ub->nnz)
  {
    bool valid_c  = (k_c < c->nnz);
    bool valid_lb = (k_lb < lb->nnz);
    bool valid_ub = (k_ub < ub->nnz);

    int i_c  = valid_c ? c->indices[k_c] : num_constraints + 1;
    int i_lb = valid_lb ? lb->indices[k_lb] : num_constraints + 1;
    int i_ub = valid_ub ? ub->indices[k_ub] : num_constraints + 1;

    int i_combined;

    i_combined = SLEQP_MIN(i_lb, i_ub);
    i_combined = SLEQP_MIN(i_combined, i_c);

    valid_c  = valid_c && (i_c == i_combined);
    valid_lb = valid_lb && (i_lb == i_combined);
    valid_ub = valid_ub && (i_ub == i_combined);

    double c_val  = valid_c ? c->data[k_c] : 0.;
    double lb_val = valid_lb ? lb->data[k_lb] : 0.;
    double ub_val = valid_ub ? ub->data[k_ub] : 0.;

    c_val = SLEQP_MAX(SLEQP_MIN(c_val, ub_val), lb_val);

    if (c_val != 0.)
    {
      SLEQP_CALL(sleqp_vec_push(result, num_variables + i_combined, c_val));
    }

    if (valid_c)
    {
      ++k_c;
    }

    if (valid_lb)
    {
      ++k_lb;
    }

    if (valid_ub)
    {
      ++k_ub;
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_restoration_problem_restore(SleqpProblem* problem,
                                  const SleqpVec* input,
                                  SleqpVec* result)
{
  const int num_variables   = sleqp_problem_num_vars(problem);
  const int num_constraints = sleqp_problem_num_cons(problem);

  assert(input->dim == (num_variables + num_constraints));
  assert(result->dim == num_variables);

  SLEQP_CALL(sleqp_vec_clear(result));

  SLEQP_CALL(sleqp_vec_reserve(result, SLEQP_MIN(input->dim, num_variables)));

  for (int k = 0; k < input->nnz; ++k)
  {
    const int j = input->indices[k];

    if (j >= num_variables)
    {
      break;
    }

    SLEQP_CALL(sleqp_vec_push(result, j, input->data[k]));
  }

  return SLEQP_OKAY;
}
