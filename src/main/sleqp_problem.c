#include "sleqp_problem.h"

#include <assert.h>
#include <math.h>

#include "sleqp_assert.h"
#include "sleqp_cmp.h"
#include "sleqp_mem.h"

typedef struct SleqpProblem
{
  int refcount;

  SleqpFunc* func;
  SleqpParams* params;

  SleqpSparseVec* var_lb;
  SleqpSparseVec* var_ub;

  SleqpSparseVec* cons_lb;
  SleqpSparseVec* cons_ub;

  SleqpSparseVec* general_lb;
  SleqpSparseVec* general_ub;

  SleqpSparseMatrix* linear_coeffs;

  SleqpSparseVec* linear_lb;
  SleqpSparseVec* linear_ub;

  int num_variables;
  int num_constraints;

  int num_linear_constraints;

  int num_general_constraints;

  SleqpSparseVec* primal;
  double* dense_cache;

  SleqpSparseVec* general_cons_val;
  SleqpSparseVec* linear_cons_val;

  SleqpSparseMatrix* general_cons_jac;

} SleqpProblem;

static SLEQP_RETCODE map_pos_inf(SleqpSparseVec* vec, double value)
{
  for(int k = 0; k < vec->nnz; ++k)
  {
    if(isinf(vec->data[k]))
    {
      vec->data[k] = value;
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE map_neg_inf(SleqpSparseVec* vec, double value)
{
  for(int k = 0; k < vec->nnz; ++k)
  {
    if(isinf(-vec->data[k]))
    {
      vec->data[k] = value;
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE check_bounds(const SleqpSparseVec* lb,
                                  const SleqpSparseVec* ub,
                                  bool cons_bounds)
{
  assert(lb->dim == ub->dim);

  int k_lb = 0, k_ub = 0;

  while(k_lb < lb->nnz || k_ub < ub->nnz)
  {
    double lb_val = 0;
    double ub_val = 0;

    bool valid_lb = (k_lb < lb->nnz);
    bool valid_ub = (k_ub < ub->nnz);

    int lb_i = valid_lb ? lb->indices[k_lb] : lb->dim + 1;
    int ub_i = valid_ub ? ub->indices[k_ub] : ub->dim + 1;

    int idx = SLEQP_MIN(lb_i, ub_i);

    valid_lb = valid_lb && (lb_i == idx);
    valid_ub = valid_ub && (ub_i == idx);

    if(valid_lb)
    {
      lb_val = lb->data[k_lb];
    }

    if(valid_ub)
    {
      ub_val = ub->data[k_ub];
    }

    if(lb_val > ub_val)
    {
      if(cons_bounds)
      {
        sleqp_log_error("Inconsistent constraint bound values at index %d: lower = %f > %f = upper",
                        idx,
                        lb_val,
                        ub_val);
      }
      else
      {
        sleqp_log_error("Inconsistent variable bound values at index %d: lower = %f > %f = upper",
                        idx,
                        lb_val,
                        ub_val);
      }

      return SLEQP_ILLEGAL_ARGUMENT;
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

static SLEQP_RETCODE problem_create(SleqpProblem** star,
                                    SleqpFunc* func,
                                    SleqpParams* params,
                                    const SleqpSparseVec* var_lb,
                                    const SleqpSparseVec* var_ub,
                                    const SleqpSparseVec* cons_lb,
                                    const SleqpSparseVec* cons_ub)
{
  assert(var_lb->dim == var_ub->dim);
  assert(cons_lb->dim == cons_ub->dim);

  sleqp_assert_msg(sleqp_sparse_vector_is_valid(var_lb), "Invalid variable bounds");
  sleqp_assert_msg(sleqp_sparse_vector_is_valid(var_ub), "Invalid variable bounds");

  sleqp_assert_msg(sleqp_sparse_vector_is_valid(cons_lb), "Invalid constraint bounds");
  sleqp_assert_msg(sleqp_sparse_vector_is_valid(cons_ub), "Invalid constraint bounds");

  const int num_constraints = sleqp_func_get_num_constraints(func);
  const int num_variables = sleqp_func_get_num_variables(func);

  SLEQP_CALL(sleqp_malloc(star));

  SleqpProblem* problem = *star;

  *problem = (SleqpProblem) {0};

  problem->refcount = 1;

  problem->num_variables = num_variables;
  problem->num_constraints = num_constraints;
  problem->num_general_constraints = num_constraints;
  problem->num_linear_constraints = 0;

  sleqp_assert_msg(var_lb->dim == num_variables, "Inconsistent variable dimensions");
  sleqp_assert_msg(var_ub->dim == num_variables, "Inconsistent variable dimensions");

  sleqp_assert_msg(cons_lb->dim == num_constraints, "Inconsistent constraint dimensions");
  sleqp_assert_msg(cons_ub->dim == num_constraints, "Inconsistent constraint dimensions");

  SLEQP_CALL(sleqp_func_capture(func));

  problem->func = func;

  SLEQP_CALL(sleqp_params_capture(params));
  problem->params = params;

  SLEQP_CALL(sleqp_sparse_vector_create(&problem->var_lb,
                                        problem->num_variables,
                                        var_lb->nnz));

  SLEQP_CALL(sleqp_sparse_vector_create(&problem->var_ub,
                                        problem->num_variables,
                                        var_ub->nnz));

  SLEQP_CALL(sleqp_sparse_vector_create(&problem->general_lb,
                                        problem->num_constraints,
                                        cons_lb->nnz));

  SLEQP_CALL(sleqp_sparse_vector_create(&problem->general_ub,
                                        problem->num_constraints,
                                        cons_ub->nnz));

  SLEQP_CALL(sleqp_sparse_vector_copy(var_lb, problem->var_lb));
  SLEQP_CALL(sleqp_sparse_vector_copy(var_ub, problem->var_ub));

  SLEQP_CALL(sleqp_sparse_vector_copy(cons_lb, problem->general_lb));
  SLEQP_CALL(sleqp_sparse_vector_copy(cons_ub, problem->general_ub));

  const double inf = sleqp_infinity();

  SLEQP_CALL(map_pos_inf(problem->var_ub, inf));
  SLEQP_CALL(map_pos_inf(problem->general_ub, inf));

  SLEQP_CALL(map_neg_inf(problem->var_lb, -inf));
  SLEQP_CALL(map_neg_inf(problem->general_lb, -inf));

  SLEQP_CALL(check_bounds(problem->var_lb, problem->var_ub, false));
  SLEQP_CALL(check_bounds(problem->general_lb, problem->general_ub, true));

  sleqp_assert_msg(sleqp_sparse_vector_is_finite(problem->var_lb), "Infinite variable bounds");
  sleqp_assert_msg(sleqp_sparse_vector_is_finite(problem->var_ub), "Infinite variable bounds");

  sleqp_assert_msg(sleqp_sparse_vector_is_finite(problem->general_lb), "Infinite constraint bounds");
  sleqp_assert_msg(sleqp_sparse_vector_is_finite(problem->general_ub), "Infinite constraint bounds");

  SLEQP_CALL(sleqp_sparse_matrix_create(&problem->linear_coeffs,
                                        0,
                                        problem->num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&problem->cons_lb,
                                              0));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&problem->cons_ub,
                                              0));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&problem->linear_lb,
                                              0));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&problem->linear_ub,
                                              0));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_problem_create_simple(SleqpProblem** star,
                                          SleqpFunc* func,
                                          SleqpParams* params,
                                          const SleqpSparseVec* var_lb,
                                          const SleqpSparseVec* var_ub,
                                          const SleqpSparseVec* general_lb,
                                          const SleqpSparseVec* general_ub)
{
  SLEQP_CALL(problem_create(star,
                            func,
                            params,
                            var_lb,
                            var_ub,
                            general_lb,
                            general_ub));

  SleqpProblem* problem = *star;

  SLEQP_CALL(sleqp_sparse_vector_resize(problem->cons_lb,
                                        problem->num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_copy(problem->general_lb,
                                      problem->cons_lb));

  SLEQP_CALL(sleqp_sparse_vector_resize(problem->cons_ub,
                                        problem->num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_copy(problem->general_ub,
                                      problem->cons_ub));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_problem_create(SleqpProblem** star,
                                   SleqpFunc* func,
                                   SleqpParams* params,
                                   const SleqpSparseVec* var_lb,
                                   const SleqpSparseVec* var_ub,
                                   const SleqpSparseVec* general_lb,
                                   const SleqpSparseVec* general_ub,
                                   const SleqpSparseMatrix* linear_coeffs,
                                   const SleqpSparseVec* linear_lb,
                                   const SleqpSparseVec* linear_ub)
{
  SLEQP_CALL(problem_create(star,
                            func,
                            params,
                            var_lb,
                            var_ub,
                            general_lb,
                            general_ub));

  SleqpProblem* problem = *star;

  const double inf = sleqp_infinity();

  const int num_variables = problem->num_variables;

  const int num_general_constraints = problem->num_constraints;
  const int num_linear_constraints = sleqp_sparse_matrix_get_num_rows(linear_coeffs);

  SLEQP_CALL(sleqp_sparse_matrix_resize(problem->linear_coeffs,
                                        num_linear_constraints,
                                        num_variables));

  SLEQP_CALL(sleqp_sparse_matrix_copy(linear_coeffs,
                                      problem->linear_coeffs));

  SLEQP_CALL(sleqp_sparse_vector_resize(problem->linear_lb, num_linear_constraints));

  SLEQP_CALL(sleqp_sparse_vector_copy(linear_lb,
                                      problem->linear_lb));

  SLEQP_CALL(sleqp_sparse_vector_resize(problem->linear_ub, num_linear_constraints));

  SLEQP_CALL(sleqp_sparse_vector_copy(linear_ub,
                                      problem->linear_ub));

  SLEQP_CALL(sleqp_sparse_vector_concat(general_lb, linear_lb, problem->cons_lb));
  SLEQP_CALL(sleqp_sparse_vector_concat(general_ub, linear_ub, problem->cons_ub));

  SLEQP_CALL(map_neg_inf(problem->linear_lb, -inf));
  SLEQP_CALL(map_neg_inf(problem->linear_ub, inf));

  sleqp_assert_msg(sleqp_sparse_matrix_valid(problem->linear_coeffs), "Invalid linear coefficients");

  sleqp_assert_msg(sleqp_sparse_vector_is_finite(problem->linear_lb), "Infinite linear bounds");
  sleqp_assert_msg(sleqp_sparse_vector_is_finite(problem->linear_ub), "Infinite linear bounds");

  sleqp_assert_msg(sleqp_sparse_matrix_get_num_cols(linear_coeffs) == num_variables,
                   "Inconsistent constraint dimensions");

  problem->num_constraints = num_general_constraints + num_linear_constraints;
  problem->num_linear_constraints = num_linear_constraints;
  problem->num_general_constraints = num_general_constraints;

  sleqp_assert_msg(linear_lb->dim == num_linear_constraints, "Inconsistent linear constraint dimensions");
  sleqp_assert_msg(linear_ub->dim == num_linear_constraints, "Inconsistent linear constraint dimensions");

  SLEQP_CALL(sleqp_sparse_vector_create_full(&problem->primal,
                                             num_variables));

  SLEQP_CALL(sleqp_alloc_array(&problem->dense_cache, num_linear_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_full(&problem->general_cons_val,
                                             num_general_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_full(&problem->linear_cons_val,
                                             num_linear_constraints));

  SLEQP_CALL(sleqp_sparse_matrix_create(&problem->general_cons_jac,
                                        num_general_constraints,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_concat(general_lb,
                                        linear_lb,
                                        problem->cons_lb));

  SLEQP_CALL(sleqp_sparse_vector_concat(general_ub,
                                        linear_ub,
                                        problem->cons_ub));

  SLEQP_CALL(check_bounds(problem->cons_lb, problem->cons_ub, true));

  return SLEQP_OKAY;
}

int sleqp_problem_num_variables(SleqpProblem* problem)
{
  return problem->num_variables;
}

SleqpSparseVec* sleqp_problem_var_lb(SleqpProblem* problem)
{
  return problem->var_lb;
}

SleqpSparseVec* sleqp_problem_var_ub(SleqpProblem* problem)
{
  return problem->var_ub;
}

SleqpSparseVec* sleqp_problem_general_lb(SleqpProblem* problem)
{
  return problem->general_lb;
}

SleqpSparseVec* sleqp_problem_general_ub(SleqpProblem* problem)
{
  return problem->general_ub;
}

SleqpSparseMatrix* sleqp_problem_linear_coeffs(SleqpProblem* problem)
{
  return problem->linear_coeffs;
}

SleqpSparseVec* sleqp_problem_linear_lb(SleqpProblem* problem)
{
  return problem->linear_lb;
}

SleqpSparseVec* sleqp_problem_linear_ub(SleqpProblem* problem)
{
  return problem->linear_ub;
}

SleqpSparseVec* sleqp_problem_cons_lb(SleqpProblem* problem)
{
  return problem->cons_lb;
}

SleqpSparseVec* sleqp_problem_cons_ub(SleqpProblem* problem)
{
  return problem->cons_ub;
}

int sleqp_problem_num_constraints(SleqpProblem* problem)
{
  return problem->num_constraints;
}

int sleqp_problem_num_linear_constraints(SleqpProblem* problem)
{
  return problem->num_linear_constraints;
}

int sleqp_problem_num_general_constraints(SleqpProblem* problem)
{
  return problem->num_general_constraints;
}

SleqpFunc* sleqp_problem_func(SleqpProblem* problem)
{
  return problem->func;
}

SLEQP_RETCODE sleqp_problem_set_value(SleqpProblem* problem,
                                      SleqpSparseVec* x,
                                      SLEQP_VALUE_REASON reason,
                                      int* func_grad_nnz,
                                      int* cons_val_nnz,
                                      int* cons_jac_nnz)
{
  SLEQP_CALL(sleqp_func_set_value(problem->func,
                                  x,
                                  reason,
                                  func_grad_nnz,
                                  cons_val_nnz,
                                  cons_jac_nnz));

  if(problem->primal)
  {
    SLEQP_CALL(sleqp_sparse_vector_copy(x, problem->primal));
  }

  if(problem->general_cons_jac)
  {
    SLEQP_CALL(sleqp_sparse_matrix_reserve(problem->general_cons_jac,
                                           *cons_jac_nnz));
  }

  (*cons_jac_nnz) += sleqp_sparse_matrix_get_nnz(problem->linear_coeffs);

  (*cons_val_nnz) += problem->num_linear_constraints;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_problem_eval(SleqpProblem* problem,
                                 const SleqpSparseVec* cons_indices,
                                 double* func_val,
                                 SleqpSparseVec* func_grad,
                                 SleqpSparseVec* cons_val,
                                 SleqpSparseMatrix* cons_jac)
{
  if(func_val)
  {
    SLEQP_CALL(sleqp_problem_val(problem, func_val));
  }

  if(func_grad)
  {
    SLEQP_CALL(sleqp_problem_grad(problem, func_grad));
  }

  if(cons_val)
  {
    SLEQP_CALL(sleqp_problem_cons_val(problem, cons_indices, cons_val));
  }

  if(cons_jac)
  {
    SLEQP_CALL(sleqp_problem_cons_jac(problem, cons_indices, cons_jac));
  }

  if(func_grad && !sleqp_sparse_vector_is_valid(func_grad))
  {
    sleqp_log_error("Function returned invalid gradient");
    return SLEQP_ILLEGAL_ARGUMENT;
  }

  if(cons_val && !sleqp_sparse_vector_is_valid(cons_val))
  {
    sleqp_log_error("Function returned invalid constraint values");
    return SLEQP_ILLEGAL_ARGUMENT;
  }

  if(cons_jac && !sleqp_sparse_matrix_valid(cons_jac))
  {
    sleqp_log_error("Function returned invalid constraint Jacobian");
    return SLEQP_ILLEGAL_ARGUMENT;
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_problem_val(SleqpProblem* problem,
                                double* func_val)
{
  return sleqp_func_val(problem->func, func_val);
}

SLEQP_RETCODE sleqp_problem_grad(SleqpProblem* problem,
                                 SleqpSparseVec* func_grad)
{
  return sleqp_func_grad(problem->func, func_grad);
}

SLEQP_RETCODE sleqp_problem_cons_val(SleqpProblem* problem,
                                     const SleqpSparseVec* cons_indices,
                                     SleqpSparseVec* cons_val)
{
  const double zero_eps = sleqp_params_get(problem->params,
                                           SLEQP_PARAM_ZERO_EPS);

  if(problem->num_linear_constraints == 0)
  {
    return sleqp_func_cons_val(problem->func,
                               cons_indices,
                               cons_val);
  }

  SLEQP_CALL(sleqp_sparse_matrix_vector_product(problem->linear_coeffs,
                                                problem->primal,
                                                problem->dense_cache));

  if(problem->num_general_constraints == 0)
  {
    return sleqp_sparse_vector_from_raw(cons_val,
                                        problem->dense_cache,
                                        problem->num_constraints,
                                        zero_eps);
  }
  else
  {
    SLEQP_CALL(sleqp_sparse_vector_from_raw(problem->linear_cons_val,
                                            problem->dense_cache,
                                            problem->num_linear_constraints,
                                            zero_eps));

    SLEQP_CALL(sleqp_func_cons_val(problem->func,
                                   cons_indices,
                                   problem->general_cons_val));

    return sleqp_sparse_vector_concat(problem->general_cons_val,
                                      problem->linear_cons_val,
                                      cons_val);
  }
}

SLEQP_RETCODE sleqp_problem_cons_jac(SleqpProblem* problem,
                                     const SleqpSparseVec* cons_indices,
                                     SleqpSparseMatrix* cons_jac)
{
  if(problem->num_linear_constraints == 0)
  {
    return sleqp_func_cons_jac(problem->func,
                               cons_indices,
                               cons_jac);
  }

  if(problem->num_general_constraints == 0)
  {
    return sleqp_sparse_matrix_copy(problem->linear_coeffs, cons_jac);
  }

  SLEQP_CALL(sleqp_func_cons_jac(problem->func,
                                 cons_indices,
                                 problem->general_cons_jac));

  return sleqp_sparse_matrix_vstack(problem->general_cons_jac,
                                    problem->linear_coeffs,
                                    cons_jac);
}

SLEQP_RETCODE sleqp_problem_hess_prod(SleqpProblem* problem,
                                      const double* func_dual,
                                      const SleqpSparseVec* direction,
                                      const SleqpSparseVec* cons_duals,
                                      SleqpSparseVec* product)
{
  return sleqp_func_hess_prod(problem->func,
                              func_dual,
                              direction,
                              cons_duals,
                              product);
}

SLEQP_RETCODE sleqp_problem_hess_bilinear(SleqpProblem* problem,
                                          const double* func_dual,
                                          const SleqpSparseVec* direction,
                                          const SleqpSparseVec* cons_duals,
                                          double* bilinear_prod)
{
  return sleqp_func_hess_bilinear(problem->func,
                                  func_dual,
                                  direction,
                                  cons_duals,
                                  bilinear_prod);
}

static SLEQP_RETCODE problem_free(SleqpProblem** star)
{
  SleqpProblem* problem = *star;

  if(!problem)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_sparse_matrix_release(&problem->general_cons_jac));

  SLEQP_CALL(sleqp_sparse_vector_free(&problem->general_cons_val));

  SLEQP_CALL(sleqp_sparse_vector_free(&problem->linear_cons_val));

  sleqp_free(&problem->dense_cache);

  SLEQP_CALL(sleqp_sparse_vector_free(&problem->primal));


  SLEQP_CALL(sleqp_sparse_vector_free(&problem->linear_ub));
  SLEQP_CALL(sleqp_sparse_vector_free(&problem->linear_lb));

  SLEQP_CALL(sleqp_sparse_matrix_release(&problem->linear_coeffs));

  SLEQP_CALL(sleqp_sparse_vector_free(&problem->cons_ub));
  SLEQP_CALL(sleqp_sparse_vector_free(&problem->cons_lb));

  SLEQP_CALL(sleqp_sparse_vector_free(&problem->var_ub));
  SLEQP_CALL(sleqp_sparse_vector_free(&problem->var_lb));

  SLEQP_CALL(sleqp_params_release(&problem->params));

  SLEQP_CALL(sleqp_func_release(&problem->func));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_problem_capture(SleqpProblem* problem)
{
  ++problem->refcount;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_problem_release(SleqpProblem** star)
{
  SleqpProblem* problem = *star;

  if(!problem)
  {
    return SLEQP_OKAY;
  }

  if(--problem->refcount == 0)
  {
    SLEQP_CALL(problem_free(star));
  }

  *star = NULL;
  
  return SLEQP_OKAY;
}
