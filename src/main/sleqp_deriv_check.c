#include "sleqp_deriv_check.h"

#include "sleqp_cmp.h"
#include "sleqp_mem.h"

struct SleqpDerivCheckData
{
  SleqpProblem* problem;
  SleqpParams* params;

  SleqpSparseVec* value_diff;

  SleqpSparseVec* hessian_estimate;

  SleqpSparseVec* hessian_left;
  SleqpSparseVec* hessian_right;

  SleqpSparseVec* hessian_prod;

  SleqpSparseVec* hessian_prod_zero;
  SleqpSparseVec* hessian_prod_cache;

  SleqpSparseVec* cons_grad_iterate;
  SleqpSparseVec* cons_grad_check_iterate;

  SleqpSparseVec* multipliers;
  SleqpSparseVec* multipliers_zero;

  SleqpIterate* iterate;
  SleqpIterate* check_iterate;
};

SLEQP_RETCODE sleqp_deriv_checker_create(SleqpDerivCheckData** star,
                                         SleqpProblem* problem,
                                         SleqpParams* params)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpDerivCheckData* data = *star;

  int num_constraints = problem->num_constraints;
  int num_variables = problem->num_variables;

  data->problem = problem;
  data->params = params;
  data->iterate = NULL;

  // this will be a unit vector
  SLEQP_CALL(sleqp_sparse_vector_create(&data->value_diff,
                                        num_variables,
                                        1));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->hessian_estimate,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->hessian_left,
                                        num_variables,
                                        1));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->hessian_right,
                                        num_variables,
                                        1));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->hessian_prod,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->hessian_prod_zero,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->hessian_prod_cache,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->cons_grad_iterate,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->cons_grad_check_iterate,
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->multipliers,
                                        num_constraints,
                                        1));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->multipliers_zero,
                                        num_constraints,
                                        1));

  SLEQP_CALL(sleqp_iterate_create(&data->check_iterate,
                                  problem,
                                  problem->var_lb));


  return SLEQP_OKAY;
}

static double get_perturbation(double perturbation_param,
                               SleqpSparseVec* x,
                               int j)
{
  double* vec_ptr = sleqp_sparse_vector_at(x, j);

  double value = vec_ptr ? (*vec_ptr) : 0.;

  value = SLEQP_ABS(value);

  value = SLEQP_MAX(value, 1.);

  return value * perturbation_param;
}

static SLEQP_RETCODE check_func_first_order_at(SleqpDerivCheckData* data,
                                               SleqpIterate* iterate,
                                               int j,
                                               bool* valid)
{
  SleqpIterate* check_iterate = data->check_iterate;
  SleqpFunc* func = data->problem->func;

  int func_grad_nnz = 0;
  int cons_val_nnz = 0;
  int cons_jac_nnz = 0;

  SleqpSparseVec* value_diff = data->value_diff;

  const double tolerance = sleqp_params_get_deriv_tolerance(data->params);

  const double pertubation = get_perturbation(sleqp_params_get_deriv_pertubation(data->params),
                                              iterate->x,
                                              j);

  SLEQP_CALL(sleqp_sparse_vector_clear(value_diff));

  SLEQP_CALL(sleqp_sparse_vector_push(value_diff, j, 1.));

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(iterate->x,
                                            value_diff,
                                            1.,
                                            pertubation,
                                            0.,
                                            check_iterate->x));

  SLEQP_CALL(sleqp_func_set_value(func,
                                  check_iterate->x,
                                  &func_grad_nnz,
                                  &cons_val_nnz,
                                  &cons_jac_nnz));

  SLEQP_CALL(sleqp_sparse_vector_reserve(check_iterate->func_grad, func_grad_nnz));

  SLEQP_CALL(sleqp_sparse_vector_reserve(check_iterate->cons_val, cons_val_nnz));

  SLEQP_CALL(sleqp_func_eval(func,
                             NULL,
                             &(check_iterate->func_val),
                             NULL,
                             NULL,
                             NULL));

  const double func_diff = check_iterate->func_val - iterate->func_val;

  const double actual_value = func_diff / pertubation;

  double* vec_ptr = sleqp_sparse_vector_at(iterate->func_grad, j);

  const double expected_value = vec_ptr ? (*vec_ptr) : 0.;

  if(!sleqp_eq(expected_value, actual_value, tolerance))
  {
    sleqp_log_error("Derivative check failed for objective function gradient at %d: "
                    "%.10e ~ %.10e",
                    j,
                    expected_value,
                    actual_value);

    *valid = false;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE check_cons_first_order_at(SleqpDerivCheckData* data,
                                               SleqpIterate* iterate,
                                               int j,
                                               bool* valid)
{
  SleqpIterate* check_iterate = data->check_iterate;
  SleqpProblem* problem = data->problem;
  SleqpFunc* func = problem->func;

  int func_grad_nnz = 0;
  int cons_val_nnz = 0;
  int cons_jac_nnz = 0;

  SleqpSparseVec* value_diff = data->value_diff;

  const double tolerance = sleqp_params_get_deriv_tolerance(data->params);

  const double pertubation = get_perturbation(sleqp_params_get_deriv_pertubation(data->params),
                                              iterate->x,
                                              j);

  SLEQP_CALL(sleqp_sparse_vector_clear(value_diff));

  SLEQP_CALL(sleqp_sparse_vector_push(value_diff, j, 1.));

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(iterate->x,
                                            value_diff,
                                            1.,
                                            pertubation,
                                            0.,
                                            check_iterate->x));

  SLEQP_CALL(sleqp_func_set_value(func,
                                  check_iterate->x,
                                  &func_grad_nnz,
                                  &cons_val_nnz,
                                  &cons_jac_nnz));

  SLEQP_CALL(sleqp_sparse_vector_reserve(check_iterate->cons_val, cons_val_nnz));

  SLEQP_CALL(sleqp_sparse_matrix_reserve(check_iterate->cons_jac, cons_jac_nnz));

  SLEQP_CALL(sleqp_func_eval(func,
                             NULL,
                             NULL,
                             NULL,
                             check_iterate->cons_val,
                             NULL));

  for(int i = 0; i < problem->num_constraints; ++i)
  {
    double* ptr;

    ptr = sleqp_sparse_matrix_at(iterate->cons_jac, i, j);

    const double expected_value = ptr ? (*ptr) : 0.;

    ptr = sleqp_sparse_vector_at(iterate->cons_val, i);

    const double lower_value = ptr ? (*ptr) : 0.;

    ptr = sleqp_sparse_vector_at(check_iterate->cons_val, i);

    const double upper_value = ptr ? (*ptr) : 0.;

    const double value_diff = upper_value - lower_value;

    const double actual_value = value_diff / pertubation;

    if(!sleqp_eq(expected_value, actual_value, tolerance))
    {
      sleqp_log_error("Derivative check failed for %d-th constraint gradient at index %d: "
                      "%.10e ~ %.10e",
                      i,
                      j,
                      expected_value,
                      actual_value);

      *valid = false;
    }
  }

  return SLEQP_OKAY;
}


SLEQP_RETCODE sleqp_deriv_check_first_order(SleqpDerivCheckData* data,
                                            SleqpIterate* iterate)
{
  SleqpFunc* func = data->problem->func;
  SleqpProblem* problem = data->problem;

  bool valid = true;

  for(int j = 0; j < problem->num_variables; ++j)
  {
    SLEQP_CALL(check_func_first_order_at(data, iterate, j, &valid));
    SLEQP_CALL(check_cons_first_order_at(data, iterate, j, &valid));
  }

  int func_grad_nnz = 0;
  int cons_val_nnz = 0;
  int cons_jac_nnz = 0;

  // restore original function value...
  SLEQP_CALL(sleqp_func_set_value(func,
                                  iterate->x,
                                  &func_grad_nnz,
                                  &cons_val_nnz,
                                  &cons_jac_nnz));

  if(!valid)
  {
    return SLEQP_INVALID_DERIV;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE check_func_second_order_at(SleqpDerivCheckData* data,
                                                SleqpIterate* iterate,
                                                int j,
                                                bool* valid)
{
  SleqpIterate* check_iterate = data->check_iterate;
  SleqpProblem* problem = data->problem;
  SleqpFunc* func = problem->func;

  int num_variables = problem->num_variables;

  int func_grad_nnz = 0;
  int cons_val_nnz = 0;
  int cons_jac_nnz = 0;

  SleqpSparseVec* value_diff = data->value_diff;

  const double tolerance = sleqp_params_get_deriv_tolerance(data->params);

  const double pertubation = get_perturbation(sleqp_params_get_deriv_pertubation(data->params),
                                              iterate->x,
                                              j);

  SLEQP_CALL(sleqp_sparse_vector_clear(value_diff));

  SLEQP_CALL(sleqp_sparse_vector_push(value_diff, j, 1.));

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(iterate->x,
                                            value_diff,
                                            1.,
                                            pertubation,
                                            0.,
                                            check_iterate->x));

  SLEQP_CALL(sleqp_func_set_value(func,
                                  check_iterate->x,
                                  &func_grad_nnz,
                                  &cons_val_nnz,
                                  &cons_jac_nnz));

  SLEQP_CALL(sleqp_sparse_vector_reserve(check_iterate->func_grad, func_grad_nnz));

  SLEQP_CALL(sleqp_func_eval(func,
                             NULL,
                             NULL,
                             check_iterate->func_grad,
                             NULL,
                             NULL));

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(iterate->func_grad,
                                            check_iterate->func_grad,
                                            -1.,
                                            1.,
                                            0,
                                            data->hessian_estimate));

  SLEQP_CALL(sleqp_sparse_vector_scale(data->hessian_estimate, 1./pertubation));

  SLEQP_CALL(sleqp_sparse_vector_clear(data->multipliers));

  SLEQP_CALL(sleqp_sparse_vector_clear(data->hessian_right));

  SLEQP_CALL(sleqp_sparse_vector_push(data->hessian_right, j, 1.));

  double one = 1.;

  for(int k = 0; k < num_variables; ++k)
  {
    SLEQP_CALL(sleqp_func_set_value(func,
                                    iterate->x,
                                    &func_grad_nnz,
                                    &cons_val_nnz,
                                    &cons_jac_nnz));

    SLEQP_CALL(sleqp_func_hess_product(func,
                                       &one,
                                       data->hessian_right,
                                       data->multipliers,
                                       data->hessian_prod));

    SLEQP_CALL(sleqp_sparse_vector_clear(data->hessian_left));

    SLEQP_CALL(sleqp_sparse_vector_push(data->hessian_left, k, 1.));

    double hessian_coeff;

    SLEQP_CALL(sleqp_sparse_vector_dot(data->hessian_left,
                                       data->hessian_prod,
                                       &hessian_coeff));

    double* vec_ptr = sleqp_sparse_vector_at(data->hessian_estimate, k);

    const double expected_value = hessian_coeff;

    const double actual_value = vec_ptr ? (*vec_ptr) : 0.;

    if(!sleqp_eq(expected_value, actual_value, tolerance))
    {
      sleqp_log_error("Derivative check failed for objective function hessian at (%d, %d): "
                      "%.10e ~ %.10e",
                      k,
                      j,
                      expected_value,
                      actual_value);

      *valid = false;
    }

  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE check_cons_second_order_at(SleqpDerivCheckData* data,
                                                SleqpIterate* iterate,
                                                int i,
                                                int j,
                                                bool* valid)
{
  SleqpIterate* check_iterate = data->check_iterate;
  SleqpProblem* problem = data->problem;
  SleqpFunc* func = problem->func;

  int num_variables = problem->num_variables;

  int func_grad_nnz = 0;
  int cons_val_nnz = 0;
  int cons_jac_nnz = 0;

  double one = 1.;

  SleqpSparseVec* value_diff = data->value_diff;

  const double tolerance = sleqp_params_get_deriv_tolerance(data->params);

  const double pertubation = get_perturbation(sleqp_params_get_deriv_pertubation(data->params),
                                              iterate->x,
                                              j);

  SLEQP_CALL(sleqp_sparse_vector_clear(value_diff));

  SLEQP_CALL(sleqp_sparse_vector_push(value_diff, j, 1.));

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(iterate->x,
                                            value_diff,
                                            1.,
                                            pertubation,
                                            0.,
                                            check_iterate->x));

  SLEQP_CALL(sleqp_func_set_value(func,
                                  check_iterate->x,
                                  &func_grad_nnz,
                                  &cons_val_nnz,
                                  &cons_jac_nnz));

  SLEQP_CALL(sleqp_sparse_matrix_reserve(check_iterate->cons_jac, cons_jac_nnz));

  SLEQP_CALL(sleqp_func_eval(func,
                             NULL,
                             NULL,
                             NULL,
                             NULL,
                             check_iterate->cons_jac));

  SLEQP_CALL(sleqp_sparse_vector_clear(data->multipliers));

  SLEQP_CALL(sleqp_sparse_vector_push(data->multipliers, i, 1));

  SLEQP_CALL(sleqp_sparse_matrix_trans_vector_product(iterate->cons_jac,
                                                      data->multipliers,
                                                      0.,
                                                      data->cons_grad_iterate));

  SLEQP_CALL(sleqp_sparse_matrix_trans_vector_product(check_iterate->cons_jac,
                                                      data->multipliers,
                                                      0.,
                                                      data->cons_grad_check_iterate));

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(data->cons_grad_iterate,
                                            data->cons_grad_check_iterate,
                                            -1.,
                                            1.,
                                            0,
                                            data->hessian_estimate));

  SLEQP_CALL(sleqp_sparse_vector_scale(data->hessian_estimate, 1./pertubation));

  SLEQP_CALL(sleqp_sparse_vector_clear(data->hessian_right));

  SLEQP_CALL(sleqp_sparse_vector_push(data->hessian_right, j, 1.));

  SLEQP_CALL(sleqp_func_hess_product(func,
                                     &one,
                                     data->hessian_right,
                                     data->multipliers_zero,
                                     data->hessian_prod_zero));

  for(int k = 0; k < num_variables; ++k)
  {
    SLEQP_CALL(sleqp_func_set_value(func,
                                    iterate->x,
                                    &func_grad_nnz,
                                    &cons_val_nnz,
                                    &cons_jac_nnz));

    SLEQP_CALL(sleqp_func_hess_product(func,
                                       &one,
                                       data->hessian_right,
                                       data->multipliers,
                                       data->hessian_prod_cache));

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(data->hessian_prod_cache,
                                              data->hessian_prod_zero,
                                              1.,
                                              -1.,
                                              0,
                                              data->hessian_prod));

    SLEQP_CALL(sleqp_sparse_vector_clear(data->hessian_left));

    SLEQP_CALL(sleqp_sparse_vector_push(data->hessian_left, k, 1.));

    double hessian_coeff;

    SLEQP_CALL(sleqp_sparse_vector_dot(data->hessian_left,
                                       data->hessian_prod,
                                       &hessian_coeff));

    double* vec_ptr = sleqp_sparse_vector_at(data->hessian_estimate, k);

    const double expected_value = hessian_coeff;

    const double actual_value = vec_ptr ? (*vec_ptr) : 0.;

    if(!sleqp_eq(expected_value, actual_value, tolerance))
    {
      sleqp_log_error("Derivative check failed for %d-th constraint hessian at (%d, %d): "
                      "%.10e ~ %.10e",
                      i,
                      k,
                      j,
                      expected_value,
                      actual_value);

      *valid = false;
    }

  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_deriv_check_second_order(SleqpDerivCheckData* data,
                                             SleqpIterate* iterate)
{
  SleqpFunc* func = data->problem->func;
  SleqpProblem* problem = data->problem;

  bool valid = true;

  for(int j = 0; j < problem->num_variables; ++j)
  {
    SLEQP_CALL(check_func_second_order_at(data, iterate, j, &valid));
  }

  for(int i = 0; i < problem->num_constraints; ++i)
  {
    for(int j = 0; j < problem->num_variables; ++j)
    {
      SLEQP_CALL(check_cons_second_order_at(data, iterate, i, j, &valid));
    }
  }

  int func_grad_nnz = 0;
  int cons_val_nnz = 0;
  int cons_jac_nnz = 0;

  // restore original function value...
  SLEQP_CALL(sleqp_func_set_value(func,
                                  iterate->x,
                                  &func_grad_nnz,
                                  &cons_val_nnz,
                                  &cons_jac_nnz));

  if(!valid)
  {
    return SLEQP_INVALID_DERIV;
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_deriv_checker_free(SleqpDerivCheckData** star)
{
  SleqpDerivCheckData* data = *star;

  if(!data)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_iterate_free(&data->check_iterate));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->multipliers_zero));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->multipliers));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->cons_grad_check_iterate));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->cons_grad_iterate));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->hessian_prod_cache));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->hessian_prod_zero));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->hessian_prod));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->hessian_right));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->hessian_left));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->hessian_estimate));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->value_diff));

  sleqp_free(star);

  return SLEQP_OKAY;
}
