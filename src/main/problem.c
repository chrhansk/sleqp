#include "problem.h"

#include <assert.h>
#include <math.h>

#include "cmp.h"
#include "error.h"
#include "fail.h"
#include "mem.h"

#include "sparse/sparse_matrix.h"

typedef struct SleqpProblem
{
  int refcount;

  SleqpFunc* func;
  SleqpParams* params;

  SleqpVec* var_lb;
  SleqpVec* var_ub;

  SleqpVec* cons_lb;
  SleqpVec* cons_ub;

  SleqpVec* general_lb;
  SleqpVec* general_ub;

  SleqpSparseMatrix* linear_coeffs;

  SleqpVec* linear_lb;
  SleqpVec* linear_ub;

  int num_variables;

  int num_linear_constraints;
  int num_general_constraints;

  SleqpVec* primal;
  double* dense_cache;

  SleqpVec* general_cons_val;
  SleqpVec* linear_cons_val;

  SleqpVec* general_cons_duals;

  SleqpSparseMatrix* general_cons_jac;

} SleqpProblem;

static SLEQP_RETCODE
check_bounds(const SleqpVec* lb, const SleqpVec* ub, bool cons_bounds)
{
  assert(lb->dim == ub->dim);

  int k_lb = 0, k_ub = 0;

  while (k_lb < lb->nnz || k_ub < ub->nnz)
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

    if (valid_lb)
    {
      lb_val = lb->data[k_lb];
    }

    if (valid_ub)
    {
      ub_val = ub->data[k_ub];
    }

    if (lb_val > ub_val)
    {
      if (cons_bounds)
      {
        sleqp_raise(SLEQP_ILLEGAL_ARGUMENT,
                    "Inconsistent constraint bound values at index %d: "
                    "lower = %f > %f = upper",
                    idx,
                    lb_val,
                    ub_val);
      }
      else
      {
        sleqp_raise(SLEQP_ILLEGAL_ARGUMENT,
                    "Inconsistent variable bound values at index %d: lower "
                    "= %f > %f = upper",
                    idx,
                    lb_val,
                    ub_val);
      }
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

static SLEQP_RETCODE
vector_map(SleqpVec* vector, double from, double to)
{
  for (int k = 0; k < vector->nnz; ++k)
  {
    if (vector->data[k] == from)
    {
      vector->data[k] = to;
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
copy_create_vector(const SleqpVec* source,
                   int expected_dim,
                   SleqpVec** target_star)
{
  sleqp_assert_msg(sleqp_vec_is_valid(source), "Invalid lower bounds");

  sleqp_assert_msg(source->dim == expected_dim, "Inconsistent dimensions");

  if (*target_star)
  {
    SLEQP_CALL(sleqp_vec_resize(*target_star, expected_dim));

    SLEQP_CALL(sleqp_vec_reserve(*target_star, source->nnz));
  }
  else
  {
    SLEQP_CALL(sleqp_vec_create(target_star, expected_dim, source->nnz));
  }

  SleqpVec* target = *target_star;

  SLEQP_CALL(sleqp_vec_copy(source, target));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
convert_lb(const SleqpVec* source_lb,
           int expected_dim,
           SleqpVec** target_lb_star)
{
  const double inf = sleqp_infinity();

  SLEQP_CALL(copy_create_vector(source_lb, expected_dim, target_lb_star));

  SleqpVec* target_lb = *target_lb_star;

  SLEQP_CALL(vector_map(target_lb, -INFINITY, -inf));

  sleqp_assert_msg(sleqp_vec_is_finite(target_lb), "Lower bound is infinite");

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
convert_ub(const SleqpVec* source_ub,
           int expected_dim,
           SleqpVec** target_ub_star)
{
  const double inf = sleqp_infinity();

  SLEQP_CALL(copy_create_vector(source_ub, expected_dim, target_ub_star));

  SleqpVec* target_ub = *target_ub_star;

  SLEQP_CALL(vector_map(target_ub, INFINITY, inf));

  sleqp_assert_msg(sleqp_vec_is_finite(target_ub), "Upper bound is infinite");

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
stack_bounds(SleqpProblem* problem)
{
  SLEQP_CALL(sleqp_vec_concat(problem->general_lb,
                              problem->linear_lb,
                              problem->cons_lb));

  SLEQP_CALL(sleqp_vec_concat(problem->general_ub,
                              problem->linear_ub,
                              problem->cons_ub));

  SLEQP_CALL(check_bounds(problem->cons_lb, problem->cons_ub, true));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
problem_create(SleqpProblem** star,
               SleqpFunc* func,
               SleqpParams* params,
               const SleqpVec* var_lb,
               const SleqpVec* var_ub,
               const SleqpVec* general_lb,
               const SleqpVec* general_ub)
{
  const int num_general   = sleqp_func_num_cons(func);
  const int num_variables = sleqp_func_num_vars(func);

  SLEQP_CALL(sleqp_malloc(star));

  SleqpProblem* problem = *star;

  *problem = (SleqpProblem){0};

  problem->refcount = 1;

  problem->num_variables           = num_variables;
  problem->num_general_constraints = num_general;

  SLEQP_CALL(convert_lb(var_lb, num_variables, &problem->var_lb));
  SLEQP_CALL(convert_ub(var_ub, num_variables, &problem->var_ub));

  SLEQP_CALL(convert_lb(general_lb, num_general, &problem->general_lb));
  SLEQP_CALL(convert_ub(general_ub, num_general, &problem->general_ub));

  problem->func = func;
  SLEQP_CALL(sleqp_func_capture(problem->func));

  problem->params = params;
  SLEQP_CALL(sleqp_params_capture(problem->params));

  SLEQP_CALL(sleqp_sparse_matrix_create(&problem->linear_coeffs,
                                        0,
                                        problem->num_variables,
                                        0));

  SLEQP_CALL(sleqp_vec_create_empty(&problem->cons_lb, 0));

  SLEQP_CALL(sleqp_vec_create_empty(&problem->cons_ub, 0));

  SLEQP_CALL(sleqp_vec_create_empty(&problem->linear_lb, 0));

  SLEQP_CALL(sleqp_vec_create_empty(&problem->linear_ub, 0));

  SLEQP_CALL(check_bounds(problem->var_lb, problem->var_ub, false));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_problem_create_simple(SleqpProblem** star,
                            SleqpFunc* func,
                            SleqpParams* params,
                            const SleqpVec* var_lb,
                            const SleqpVec* var_ub,
                            const SleqpVec* general_lb,
                            const SleqpVec* general_ub)
{
  SLEQP_CALL(
    problem_create(star, func, params, var_lb, var_ub, general_lb, general_ub));

  SleqpProblem* problem = *star;

  SLEQP_CALL(stack_bounds(problem));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_problem_create(SleqpProblem** star,
                     SleqpFunc* func,
                     SleqpParams* params,
                     const SleqpVec* var_lb,
                     const SleqpVec* var_ub,
                     const SleqpVec* general_lb,
                     const SleqpVec* general_ub,
                     const SleqpSparseMatrix* linear_coeffs,
                     const SleqpVec* linear_lb,
                     const SleqpVec* linear_ub)
{
  SLEQP_CALL(
    problem_create(star, func, params, var_lb, var_ub, general_lb, general_ub));

  SleqpProblem* problem = *star;

  const int num_variables = problem->num_variables;

  const int num_linear_constraints
    = sleqp_sparse_matrix_num_rows(linear_coeffs);

  {
    SLEQP_CALL(sleqp_sparse_matrix_resize(problem->linear_coeffs,
                                          num_linear_constraints,
                                          num_variables));

    SLEQP_CALL(sleqp_sparse_matrix_copy(linear_coeffs, problem->linear_coeffs));

    sleqp_assert_msg(sleqp_sparse_matrix_is_valid(problem->linear_coeffs),
                     "Linear coefficient matrix is invalid");

    sleqp_assert_msg(sleqp_sparse_matrix_is_finite(problem->linear_coeffs),
                     "Linear coefficient matrix is not all-finite");

    sleqp_assert_msg(sleqp_sparse_matrix_num_cols(problem->linear_coeffs)
                       == num_variables,
                     "Linear constraint dimensions are inconsistent");
  }

  SLEQP_CALL(
    convert_lb(linear_lb, num_linear_constraints, &problem->linear_lb));
  SLEQP_CALL(
    convert_ub(linear_ub, num_linear_constraints, &problem->linear_ub));

  problem->num_linear_constraints = num_linear_constraints;

  {
    SLEQP_CALL(sleqp_sparse_matrix_create(&problem->general_cons_jac,
                                          problem->num_general_constraints,
                                          num_variables,
                                          0));

    SLEQP_CALL(sleqp_vec_create_full(&problem->primal, num_variables));

    SLEQP_CALL(sleqp_vec_create_empty(&problem->general_cons_val,
                                      problem->num_general_constraints));

    SLEQP_CALL(sleqp_vec_create_empty(&problem->linear_cons_val,
                                      num_linear_constraints));

    SLEQP_CALL(sleqp_vec_create_empty(&problem->general_cons_duals,
                                      problem->num_general_constraints));

    SLEQP_CALL(
      sleqp_alloc_array(&problem->dense_cache, num_linear_constraints));
  }

  SLEQP_CALL(stack_bounds(problem));

  return SLEQP_OKAY;
}

bool
sleqp_problem_has_nonlinear_cons(SleqpProblem* problem)
{
  return problem->num_general_constraints != 0;
}

int
sleqp_problem_num_vars(SleqpProblem* problem)
{
  return problem->num_variables;
}

SleqpVec*
sleqp_problem_vars_lb(SleqpProblem* problem)
{
  return problem->var_lb;
}

SleqpVec*
sleqp_problem_vars_ub(SleqpProblem* problem)
{
  return problem->var_ub;
}

SleqpVec*
sleqp_problem_general_lb(SleqpProblem* problem)
{
  return problem->general_lb;
}

SleqpVec*
sleqp_problem_general_ub(SleqpProblem* problem)
{
  return problem->general_ub;
}

SleqpSparseMatrix*
sleqp_problem_linear_coeffs(SleqpProblem* problem)
{
  return problem->linear_coeffs;
}

SleqpVec*
sleqp_problem_linear_lb(SleqpProblem* problem)
{
  return problem->linear_lb;
}

SleqpVec*
sleqp_problem_linear_ub(SleqpProblem* problem)
{
  return problem->linear_ub;
}

SleqpVec*
sleqp_problem_cons_lb(SleqpProblem* problem)
{
  return problem->cons_lb;
}

SleqpVec*
sleqp_problem_cons_ub(SleqpProblem* problem)
{
  return problem->cons_ub;
}

int
sleqp_problem_num_cons(SleqpProblem* problem)
{
  return problem->num_general_constraints + problem->num_linear_constraints;
}

int
sleqp_problem_num_lin_cons(SleqpProblem* problem)
{
  return problem->num_linear_constraints;
}

int
sleqp_problem_num_gen_cons(SleqpProblem* problem)
{
  return problem->num_general_constraints;
}

SleqpFunc*
sleqp_problem_func(SleqpProblem* problem)
{
  return problem->func;
}

SLEQP_RETCODE
sleqp_problem_set_value(SleqpProblem* problem,
                        SleqpVec* x,
                        SLEQP_VALUE_REASON reason,
                        bool* reject,
                        int* obj_grad_nnz,
                        int* cons_val_nnz,
                        int* cons_jac_nnz)
{
  SLEQP_CALL(sleqp_func_set_value(problem->func,
                                  x,
                                  reason,
                                  reject,
                                  obj_grad_nnz,
                                  cons_val_nnz,
                                  cons_jac_nnz));

  if (problem->primal)
  {
    SLEQP_CALL(sleqp_vec_copy(x, problem->primal));
  }

  if (problem->general_cons_jac)
  {
    SLEQP_CALL(
      sleqp_sparse_matrix_reserve(problem->general_cons_jac, *cons_jac_nnz));
  }

  (*cons_jac_nnz) += sleqp_sparse_matrix_nnz(problem->linear_coeffs);

  (*cons_val_nnz) += problem->num_linear_constraints;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_problem_eval(SleqpProblem* problem,
                   double* obj_val,
                   SleqpVec* obj_grad,
                   SleqpVec* cons_val,
                   SleqpSparseMatrix* cons_jac)
{
  if (obj_val)
  {
    SLEQP_CALL(sleqp_problem_obj_val(problem, obj_val));
  }

  if (obj_grad)
  {
    SLEQP_CALL(sleqp_problem_obj_grad(problem, obj_grad));
  }

  if (cons_val)
  {
    SLEQP_CALL(sleqp_problem_cons_val(problem, cons_val));
  }

  if (cons_jac)
  {
    SLEQP_CALL(sleqp_problem_cons_jac(problem, cons_jac));
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_problem_obj_val(SleqpProblem* problem, double* obj_val)
{
  return sleqp_func_obj_val(problem->func, obj_val);
}

SLEQP_RETCODE
sleqp_problem_obj_grad(SleqpProblem* problem, SleqpVec* obj_grad)
{
  return sleqp_func_obj_grad(problem->func, obj_grad);
}

SLEQP_RETCODE
sleqp_problem_cons_val(SleqpProblem* problem, SleqpVec* cons_val)
{
  const double zero_eps
    = sleqp_params_value(problem->params, SLEQP_PARAM_ZERO_EPS);

  if (problem->num_linear_constraints == 0)
  {
    return sleqp_func_cons_val(problem->func, cons_val);
  }

  SLEQP_CALL(sleqp_sparse_matrix_vector_product(problem->linear_coeffs,
                                                problem->primal,
                                                problem->dense_cache));

  if (problem->num_general_constraints == 0)
  {
    return sleqp_vec_from_raw(cons_val,
                              problem->dense_cache,
                              problem->num_linear_constraints,
                              zero_eps);
  }
  else
  {
    SLEQP_CALL(sleqp_vec_from_raw(problem->linear_cons_val,
                                  problem->dense_cache,
                                  problem->num_linear_constraints,
                                  zero_eps));

    SLEQP_CALL(sleqp_func_cons_val(problem->func, problem->general_cons_val));

    return sleqp_vec_concat(problem->general_cons_val,
                            problem->linear_cons_val,
                            cons_val);
  }
}

SLEQP_RETCODE
sleqp_problem_cons_jac(SleqpProblem* problem, SleqpSparseMatrix* cons_jac)
{
  if (problem->num_linear_constraints == 0)
  {
    return sleqp_func_cons_jac(problem->func, cons_jac);
  }

  if (problem->num_general_constraints == 0)
  {
    return sleqp_sparse_matrix_copy(problem->linear_coeffs, cons_jac);
  }

  SLEQP_CALL(sleqp_func_cons_jac(problem->func, problem->general_cons_jac));

  return sleqp_sparse_matrix_vstack(problem->general_cons_jac,
                                    problem->linear_coeffs,
                                    cons_jac);
}

static SLEQP_RETCODE
prepare_cons_duals(SleqpProblem* problem, const SleqpVec* cons_duals)
{
  SLEQP_CALL(sleqp_vec_reserve(problem->general_cons_duals, cons_duals->nnz));

  SLEQP_CALL(sleqp_vec_clear(problem->general_cons_duals));

  const int num_general = problem->num_general_constraints;

  for (int k = 0; k < cons_duals->nnz; ++k)
  {
    int i = cons_duals->indices[k];

    if (i >= num_general)
    {
      break;
    }

    SLEQP_CALL(
      sleqp_vec_push(problem->general_cons_duals, i, cons_duals->data[k]));
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_problem_hess_prod(SleqpProblem* problem,
                        const double* obj_dual,
                        const SleqpVec* direction,
                        const SleqpVec* cons_duals,
                        SleqpVec* product)
{
  if (problem->num_linear_constraints == 0)
  {
    return sleqp_func_hess_prod(problem->func,
                                obj_dual,
                                direction,
                                cons_duals,
                                product);
  }

  SLEQP_CALL(prepare_cons_duals(problem, cons_duals));

  return sleqp_func_hess_prod(problem->func,
                              obj_dual,
                              direction,
                              problem->general_cons_duals,
                              product);
}

SLEQP_RETCODE
sleqp_problem_hess_bilinear(SleqpProblem* problem,
                            const double* obj_dual,
                            const SleqpVec* direction,
                            const SleqpVec* cons_duals,
                            double* bilinear_prod)
{
  if (problem->num_linear_constraints == 0)
  {
    return sleqp_func_hess_bilinear(problem->func,
                                    obj_dual,
                                    direction,
                                    cons_duals,
                                    bilinear_prod);
  }

  SLEQP_CALL(prepare_cons_duals(problem, cons_duals));

  return sleqp_func_hess_bilinear(problem->func,
                                  obj_dual,
                                  direction,
                                  problem->general_cons_duals,
                                  bilinear_prod);
}

static bool
vec_has_finite(const SleqpVec* vec)
{
  assert(sleqp_vec_is_valid(vec));
  assert(sleqp_vec_is_finite(vec));

  if (vec->nnz != vec->dim)
  {
    return true;
  }

  for (int k = 0; k < vec->nnz; ++k)
  {
    if (sleqp_is_finite(SLEQP_ABS(vec->data[k])))
    {
      return true;
    }
  }

  return false;
}

bool
sleqp_problem_is_unconstrained(SleqpProblem* problem)
{
  return !(vec_has_finite(sleqp_problem_vars_lb(problem))
           || vec_has_finite(sleqp_problem_vars_ub(problem))
           || vec_has_finite(sleqp_problem_cons_lb(problem))
           || vec_has_finite(sleqp_problem_cons_ub(problem)));
}

static SLEQP_RETCODE
problem_free(SleqpProblem** star)
{
  SleqpProblem* problem = *star;

  if (!problem)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_sparse_matrix_release(&problem->general_cons_jac));

  SLEQP_CALL(sleqp_vec_free(&problem->general_cons_val));

  SLEQP_CALL(sleqp_vec_free(&problem->linear_cons_val));

  SLEQP_CALL(sleqp_vec_free(&problem->general_cons_duals));

  sleqp_free(&problem->dense_cache);

  SLEQP_CALL(sleqp_vec_free(&problem->primal));

  SLEQP_CALL(sleqp_vec_free(&problem->linear_ub));
  SLEQP_CALL(sleqp_vec_free(&problem->linear_lb));

  SLEQP_CALL(sleqp_sparse_matrix_release(&problem->linear_coeffs));

  SLEQP_CALL(sleqp_vec_free(&problem->general_ub));
  SLEQP_CALL(sleqp_vec_free(&problem->general_lb));

  SLEQP_CALL(sleqp_vec_free(&problem->cons_ub));
  SLEQP_CALL(sleqp_vec_free(&problem->cons_lb));

  SLEQP_CALL(sleqp_vec_free(&problem->var_ub));
  SLEQP_CALL(sleqp_vec_free(&problem->var_lb));

  SLEQP_CALL(sleqp_params_release(&problem->params));

  SLEQP_CALL(sleqp_func_release(&problem->func));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_problem_capture(SleqpProblem* problem)
{
  ++problem->refcount;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_problem_release(SleqpProblem** star)
{
  SleqpProblem* problem = *star;

  if (!problem)
  {
    return SLEQP_OKAY;
  }

  if (--problem->refcount == 0)
  {
    SLEQP_CALL(problem_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
