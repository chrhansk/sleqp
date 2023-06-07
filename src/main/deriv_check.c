#include "deriv_check.h"

#include "cmp.h"
#include "error.h"
#include "iterate.h"
#include "log.h"
#include "mem.h"
#include "problem.h"
#include "sparse/mat.h"
#include "settings.h"

struct SleqpDerivChecker
{
  SleqpProblem* problem;
  SleqpSettings* settings;

  SleqpVec* unit_direction;

  SleqpVec* hessian_estimate;

  SleqpVec* hessian_left;
  SleqpVec* hessian_right;

  SleqpVec* hessian_prod;
  SleqpVec* hessian_prod_cons;
  SleqpVec* hessian_prod_func;

  SleqpVec* hessian_prod_cache;

  SleqpMat* hessian_cons_prods;

  SleqpVec* cons_grad_iterate;
  SleqpVec* cons_grad_check_iterate;

  SleqpVec* transposed_jacobian_product;
  SleqpVec* combined_cons_grad_iterate;
  SleqpVec* combined_cons_grad_check_iterate;

  SleqpVec* multipliers;

  bool valid_deriv;

  SleqpIterate* iterate;
  SleqpIterate* check_iterate;

  SleqpVec* jac_row;
  SleqpVec* check_jac_row;
};

static SLEQP_RETCODE
restore_iterate(SleqpDerivChecker* deriv_checker, SleqpIterate* iterate)
{
  SleqpProblem* problem = deriv_checker->problem;

  bool reject = false;

  SLEQP_CALL(sleqp_problem_set_value(problem,
                                     sleqp_iterate_primal(iterate),
                                     SLEQP_VALUE_REASON_CHECKING_DERIV,
                                     &reject));

  if (reject)
  {
    sleqp_raise(SLEQP_INTERNAL_ERROR,
                "Function rejected restoration after derivative check");
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
set_check_iterate(SleqpDerivChecker* deriv_checker, bool* reject)
{
  SleqpIterate* check_iterate = deriv_checker->check_iterate;
  SleqpProblem* problem       = deriv_checker->problem;

  SLEQP_CALL(sleqp_problem_set_value(problem,
                                     sleqp_iterate_primal(check_iterate),
                                     SLEQP_VALUE_REASON_CHECKING_DERIV,
                                     reject));

  SLEQP_CALL(sleqp_iterate_reserve(check_iterate, problem));

  return SLEQP_OKAY;
}

static double
get_perturbation(const SleqpDerivChecker* deriv_checker,
                 const SleqpIterate* iterate,
                 int j)
{
  SleqpSettings* settings = deriv_checker->settings;

  SleqpVec* primal = sleqp_iterate_primal(iterate);

  double base_perturbation
    = sleqp_settings_real_value(settings, SLEQP_SETTINGS_REAL_DERIV_PERTURBATION);

  double value = sleqp_vec_value_at(primal, j);

  value = SLEQP_ABS(value);

  value = SLEQP_MAX(value, 1.);

  return value * base_perturbation;
}

static SLEQP_RETCODE
create_unit_direction(SleqpVec* direction, const int j)
{
  SLEQP_CALL(sleqp_vec_clear(direction));

  SLEQP_CALL(sleqp_vec_push(direction, j, 1.));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
create_perturbed_unit_direction(SleqpDerivChecker* deriv_checker,
                                const SleqpIterate* iterate,
                                const int j,
                                double* perturbation)
{
  SleqpVec* unit_direction = deriv_checker->unit_direction;

  *perturbation = get_perturbation(deriv_checker, iterate, j);

  SLEQP_CALL(sleqp_vec_clear(unit_direction));

  SLEQP_CALL(sleqp_vec_push(unit_direction, j, *perturbation));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_hessian_cons_products(SleqpDerivChecker* deriv_checker, int j)
{
  SleqpProblem* problem = deriv_checker->problem;

  SleqpVec* unit_direction = deriv_checker->unit_direction;

  const double zero_eps
    = sleqp_settings_real_value(deriv_checker->settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  SLEQP_CALL(create_unit_direction(unit_direction, j));

  SleqpMat* hessian_prod_matrix = deriv_checker->hessian_cons_prods;

  SLEQP_CALL(sleqp_mat_clear(hessian_prod_matrix));

  const int num_constraints = sleqp_problem_num_cons(problem);

  for (int k = 0; k < num_constraints; ++k)
  {
    SLEQP_CALL(create_unit_direction(deriv_checker->multipliers, k));

    SLEQP_CALL(sleqp_problem_hess_prod(problem,
                                       unit_direction,
                                       deriv_checker->multipliers,
                                       deriv_checker->hessian_prod_cache));

    SLEQP_CALL(sleqp_vec_add_scaled(deriv_checker->hessian_prod_cache,
                                    deriv_checker->hessian_prod_func,
                                    1.,
                                    -1.,
                                    zero_eps,
                                    deriv_checker->hessian_prod_cons));

    const int nnz = sleqp_mat_nnz(hessian_prod_matrix);

    SLEQP_CALL(sleqp_mat_reserve(hessian_prod_matrix,
                                 nnz + deriv_checker->hessian_prod_cons->nnz));

    SLEQP_CALL(sleqp_mat_push_col(hessian_prod_matrix, k));

    SLEQP_CALL(sleqp_mat_push_vec(hessian_prod_matrix,
                                  k,
                                  deriv_checker->hessian_prod_cons));
  }

  assert(sleqp_mat_is_valid(hessian_prod_matrix));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_hessian_products(SleqpDerivChecker* deriv_checker,
                         SLEQP_DERIV_CHECK flags,
                         SleqpIterate* iterate,
                         const int j)
{
  SleqpProblem* problem = deriv_checker->problem;

  SleqpVec* unit_direction = deriv_checker->unit_direction;

  SLEQP_CALL(create_unit_direction(unit_direction, j));

  if (flags & SLEQP_DERIV_CHECK_SECOND_CONS
      || flags & SLEQP_DERIV_CHECK_SECOND_OBJ)
  {
    SLEQP_CALL(sleqp_vec_clear(deriv_checker->multipliers));

    SLEQP_CALL(sleqp_problem_hess_prod(problem,
                                       unit_direction,
                                       deriv_checker->multipliers,
                                       deriv_checker->hessian_prod_func));
  }

  if (flags & SLEQP_DERIV_CHECK_SECOND_CONS)
  {
    SLEQP_CALL(compute_hessian_cons_products(deriv_checker, j));
  }

  if (flags & SLEQP_DERIV_CHECK_SECOND_SIMPLE)
  {
    SleqpVec* multipliers = sleqp_iterate_cons_dual(iterate);

    SLEQP_CALL(sleqp_problem_hess_prod(problem,
                                       unit_direction,
                                       multipliers,
                                       deriv_checker->hessian_prod));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
create_check_iterate(SleqpDerivChecker* deriv_checker,
                     const SleqpIterate* iterate,
                     const int j,
                     double* perturbation)
{
  SleqpIterate* check_iterate = deriv_checker->check_iterate;
  SleqpVec* unit_direction    = deriv_checker->unit_direction;

  const double zero_eps
    = sleqp_settings_real_value(deriv_checker->settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  SLEQP_CALL(
    create_perturbed_unit_direction(deriv_checker, iterate, j, perturbation));

  SLEQP_CALL(sleqp_vec_add(sleqp_iterate_primal(iterate),
                           unit_direction,
                           zero_eps,
                           sleqp_iterate_primal(check_iterate)));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
eval_at_check_iterate(SleqpDerivChecker* deriv_checker, SLEQP_DERIV_CHECK flags)
{
  SleqpIterate* check_iterate = deriv_checker->check_iterate;
  SleqpProblem* problem       = deriv_checker->problem;

  double obj_val;

  double* obj_ptr    = NULL;
  SleqpVec* cons_val = NULL;
  SleqpVec* obj_grad = NULL;
  SleqpMat* cons_jac = NULL;

  if (flags & SLEQP_DERIV_CHECK_FIRST_OBJ)
  {
    obj_ptr = &obj_val;
  }

  if (flags & SLEQP_DERIV_CHECK_FIRST_CONS)
  {
    cons_val = sleqp_iterate_cons_val(check_iterate);
  }

  if ((flags & SLEQP_DERIV_CHECK_SECOND_OBJ)
      || (flags & SLEQP_DERIV_CHECK_SECOND_SIMPLE))
  {
    obj_grad = sleqp_iterate_obj_grad(check_iterate);
  }

  if ((flags & SLEQP_DERIV_CHECK_SECOND_CONS)
      || (flags & SLEQP_DERIV_CHECK_SECOND_SIMPLE))
  {
    cons_jac = sleqp_iterate_cons_jac(check_iterate);
  }

  SLEQP_CALL(
    sleqp_problem_eval(problem, obj_ptr, obj_grad, cons_val, cons_jac));

  if (obj_ptr)
  {
    SLEQP_CALL(sleqp_iterate_set_obj_val(check_iterate, obj_val));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
check_deriv_func_first_order(SleqpDerivChecker* deriv_checker,
                             SleqpIterate* iterate,
                             int j,
                             double perturbation)
{
  SleqpIterate* check_iterate = deriv_checker->check_iterate;

  const double tolerance
    = sleqp_settings_real_value(deriv_checker->settings, SLEQP_SETTINGS_REAL_DERIV_TOL);

  const double func_diff
    = sleqp_iterate_obj_val(check_iterate) - sleqp_iterate_obj_val(iterate);

  const double actual_value = func_diff / perturbation;

  const double expected_value
    = sleqp_vec_value_at(sleqp_iterate_obj_grad(iterate), j);

  if (!sleqp_is_eq(expected_value, actual_value, tolerance))
  {
    sleqp_log_error(
      "Derivative check failed for objective function gradient at %d: "
      "grad = %.10e != %.10e = findiff",
      j,
      expected_value,
      actual_value);

    deriv_checker->valid_deriv = false;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
check_deriv_cons_first_order(SleqpDerivChecker* deriv_checker,
                             SleqpIterate* iterate,
                             int j,
                             double perturbation)
{
  SleqpIterate* check_iterate = deriv_checker->check_iterate;
  SleqpProblem* problem       = deriv_checker->problem;
  SleqpMat* cons_jac          = sleqp_iterate_cons_jac(iterate);

  const int num_constraints = sleqp_problem_num_cons(problem);

  const double tolerance
    = sleqp_settings_real_value(deriv_checker->settings, SLEQP_SETTINGS_REAL_DERIV_TOL);

  for (int i = 0; i < num_constraints; ++i)
  {
    const double expected_value = sleqp_mat_value_at(cons_jac, i, j);

    const double lower_value
      = sleqp_vec_value_at(sleqp_iterate_cons_val(iterate), i);

    const double upper_value
      = sleqp_vec_value_at(sleqp_iterate_cons_val(check_iterate), i);

    const double unit_direction = upper_value - lower_value;

    const double actual_value = unit_direction / perturbation;

    if (!sleqp_is_eq(expected_value, actual_value, tolerance))
    {
      sleqp_log_error(
        "Derivative check failed for %d-th constraint gradient at index %d: "
        "jac = %.10e != %.10e = findiff",
        i,
        j,
        expected_value,
        actual_value);

      deriv_checker->valid_deriv = false;
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
check_deriv_func_second_order(SleqpDerivChecker* deriv_checker,
                              SleqpIterate* iterate,
                              int j,
                              double perturbation)
{
  SleqpIterate* check_iterate = deriv_checker->check_iterate;
  SleqpProblem* problem       = deriv_checker->problem;

  const int num_variables = sleqp_problem_num_vars(problem);

  const double tolerance
    = sleqp_settings_real_value(deriv_checker->settings, SLEQP_SETTINGS_REAL_DERIV_TOL);

  SLEQP_CALL(sleqp_vec_add_scaled(sleqp_iterate_obj_grad(iterate),
                                  sleqp_iterate_obj_grad(check_iterate),
                                  -1.,
                                  1.,
                                  0,
                                  deriv_checker->hessian_estimate));

  SLEQP_CALL(
    sleqp_vec_scale(deriv_checker->hessian_estimate, 1. / perturbation));

  for (int k = 0; k < num_variables; ++k)
  {
    SLEQP_CALL(create_unit_direction(deriv_checker->hessian_left, k));

    double expected_value;

    SLEQP_CALL(sleqp_vec_dot(deriv_checker->hessian_left,
                             deriv_checker->hessian_prod_func,
                             &expected_value));

    const double actual_value
      = sleqp_vec_value_at(deriv_checker->hessian_estimate, k);

    if (!sleqp_is_eq(expected_value, actual_value, tolerance))
    {
      sleqp_log_error(
        "Derivative check failed for objective function hessian at (%d, %d): "
        "hess = %.10e != %.10e = findiff",
        k,
        j,
        expected_value,
        actual_value);

      deriv_checker->valid_deriv = false;
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
check_deriv_cons_second_order(SleqpDerivChecker* deriv_checker,
                              SleqpIterate* iterate,
                              int j,
                              double perturbation)
{
  SleqpIterate* check_iterate = deriv_checker->check_iterate;
  SleqpProblem* problem       = deriv_checker->problem;

  const int num_constraints = sleqp_problem_num_cons(problem);
  const int num_variables   = sleqp_problem_num_vars(problem);

  const double zero_eps
    = sleqp_settings_real_value(deriv_checker->settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  const double tolerance
    = sleqp_settings_real_value(deriv_checker->settings, SLEQP_SETTINGS_REAL_DERIV_TOL);

  SleqpMat* cons_jac  = sleqp_iterate_cons_jac(iterate);
  SleqpMat* check_jac = sleqp_iterate_cons_jac(check_iterate);

  SleqpVec* expected_hessian_prod = deriv_checker->hessian_prod_cons;
  SleqpVec* actual_hessian_prod   = deriv_checker->hessian_estimate;

  for (int i = 0; i < num_constraints; ++i)
  {
    SLEQP_CALL(sleqp_mat_col(deriv_checker->hessian_cons_prods,
                             i,
                             expected_hessian_prod));

    SLEQP_CALL(create_unit_direction(deriv_checker->multipliers, i));

    SLEQP_CALL(sleqp_mat_mult_vec_trans(cons_jac,
                                        deriv_checker->multipliers,
                                        zero_eps,
                                        deriv_checker->jac_row));

    SLEQP_CALL(sleqp_mat_mult_vec_trans(check_jac,
                                        deriv_checker->multipliers,
                                        zero_eps,
                                        deriv_checker->check_jac_row));

    SLEQP_CALL(sleqp_vec_add_scaled(deriv_checker->check_jac_row,
                                    deriv_checker->jac_row,
                                    1.,
                                    -1.,
                                    zero_eps,
                                    actual_hessian_prod));

    SLEQP_CALL(sleqp_vec_scale(actual_hessian_prod, 1. / perturbation));

    // validation
    for (int k = 0; k < num_variables; ++k)
    {
      const double expected_value
        = sleqp_vec_value_at(expected_hessian_prod, k);

      const double actual_value = sleqp_vec_value_at(actual_hessian_prod, k);

      if (!sleqp_is_eq(expected_value, actual_value, tolerance))
      {
        sleqp_log_error(
          "Derivative check failed for %d-th constraint hessian at (%d, %d): "
          "hess = %.10e != %.10e = findiff",
          i,
          k,
          j,
          expected_value,
          actual_value);

        deriv_checker->valid_deriv = false;
      }
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_combined_cons_grad(SleqpDerivChecker* deriv_checker,
                           SleqpVec* multipliers,
                           SleqpIterate* iterate,
                           SleqpVec* result)
{
  SleqpSettings* settings = deriv_checker->settings;

  const double zero_eps = sleqp_settings_real_value(settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  SLEQP_CALL(
    sleqp_mat_mult_vec_trans(sleqp_iterate_cons_jac(iterate),
                             multipliers,
                             zero_eps,
                             deriv_checker->transposed_jacobian_product));

  SLEQP_CALL(sleqp_vec_add(sleqp_iterate_obj_grad(iterate),
                           deriv_checker->transposed_jacobian_product,
                           zero_eps,
                           result));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
check_deriv_simple_second_order(SleqpDerivChecker* deriv_checker,
                                SleqpIterate* iterate,
                                int j,
                                double perturbation)
{
  SleqpIterate* check_iterate = deriv_checker->check_iterate;
  SleqpProblem* problem       = deriv_checker->problem;
  SleqpVec* multipliers       = sleqp_iterate_cons_dual(iterate);

  const int num_variables = sleqp_problem_num_vars(problem);

  const double tolerance
    = sleqp_settings_real_value(deriv_checker->settings, SLEQP_SETTINGS_REAL_DERIV_TOL);

  // Compute Hessian product estimate
  {
    SLEQP_CALL(
      compute_combined_cons_grad(deriv_checker,
                                 multipliers,
                                 iterate,
                                 deriv_checker->combined_cons_grad_iterate));

    SLEQP_CALL(create_check_iterate(deriv_checker, iterate, j, &perturbation));

    SLEQP_CALL(compute_combined_cons_grad(
      deriv_checker,
      multipliers,
      check_iterate,
      deriv_checker->combined_cons_grad_check_iterate));

    SLEQP_CALL(
      sleqp_vec_add_scaled(deriv_checker->combined_cons_grad_iterate,
                           deriv_checker->combined_cons_grad_check_iterate,
                           -1.,
                           1.,
                           0,
                           deriv_checker->hessian_estimate));

    SLEQP_CALL(
      sleqp_vec_scale(deriv_checker->hessian_estimate, 1. / perturbation));
  }

  // Validation
  for (int k = 0; k < num_variables; ++k)
  {
    const double actual_value
      = sleqp_vec_value_at(deriv_checker->hessian_estimate, k);

    const double expected_value
      = sleqp_vec_value_at(deriv_checker->hessian_prod, k);

    if (!sleqp_is_eq(expected_value, actual_value, tolerance))
    {
      sleqp_log_error(
        "Derivative check failed for Hessian of Lagrangean at (%d, %d): "
        "hess = %.10e != %.10e = findiff",
        k,
        j,
        expected_value,
        actual_value);

      deriv_checker->valid_deriv = false;
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
check_deriv(SleqpDerivChecker* deriv_checker,
            SLEQP_DERIV_CHECK flags,
            SleqpIterate* iterate,
            int j,
            double perturbation)
{
  if (flags & SLEQP_DERIV_CHECK_FIRST_OBJ)
  {
    SLEQP_CALL(
      check_deriv_func_first_order(deriv_checker, iterate, j, perturbation));
  }

  if (flags & SLEQP_DERIV_CHECK_FIRST_CONS)
  {
    SLEQP_CALL(
      check_deriv_cons_first_order(deriv_checker, iterate, j, perturbation));
  }

  if (flags & SLEQP_DERIV_CHECK_SECOND_OBJ)
  {
    SLEQP_CALL(
      check_deriv_func_second_order(deriv_checker, iterate, j, perturbation));
  }

  if (flags & SLEQP_DERIV_CHECK_SECOND_CONS)
  {
    SLEQP_CALL(
      check_deriv_cons_second_order(deriv_checker, iterate, j, perturbation));
  }

  if (flags & SLEQP_DERIV_CHECK_SECOND_SIMPLE)
  {
    SLEQP_CALL(
      check_deriv_simple_second_order(deriv_checker, iterate, j, perturbation));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
eval_and_check_deriv(SleqpDerivChecker* deriv_checker,
                     SleqpIterate* iterate,
                     SLEQP_DERIV_CHECK flags,
                     int j)
{
  double perturbation;

  bool reject = false;

  SLEQP_CALL(create_check_iterate(deriv_checker, iterate, j, &perturbation));

  SLEQP_CALL(set_check_iterate(deriv_checker, &reject));

  if (reject)
  {
    sleqp_log_warn("Function rejected derivative check");
  }
  else
  {
    SLEQP_CALL(compute_hessian_products(deriv_checker, flags, iterate, j));

    SLEQP_CALL(eval_at_check_iterate(deriv_checker, flags));

    SLEQP_CALL(check_deriv(deriv_checker, flags, iterate, j, perturbation));
  }

  SLEQP_CALL(restore_iterate(deriv_checker, iterate));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_deriv_checker_create(SleqpDerivChecker** deriv_checker,
                           SleqpProblem* problem,
                           SleqpSettings* settings)
{
  SLEQP_CALL(sleqp_malloc(deriv_checker));

  SleqpDerivChecker* data = *deriv_checker;

  const int num_constraints = sleqp_problem_num_cons(problem);
  const int num_variables   = sleqp_problem_num_vars(problem);

  data->problem = problem;
  SLEQP_CALL(sleqp_problem_capture(data->problem));

  SLEQP_CALL(sleqp_settings_capture(settings));
  data->settings = settings;

  data->iterate = NULL;

  // this will be a unit vector
  SLEQP_CALL(sleqp_vec_create(&data->unit_direction, num_variables, 1));

  SLEQP_CALL(sleqp_vec_create_empty(&data->hessian_estimate, num_variables));

  SLEQP_CALL(sleqp_vec_create(&data->hessian_left, num_variables, 1));

  SLEQP_CALL(sleqp_vec_create(&data->hessian_right, num_variables, 1));

  SLEQP_CALL(sleqp_vec_create_empty(&data->hessian_prod, num_variables));

  SLEQP_CALL(sleqp_vec_create_empty(&data->hessian_prod_cons, num_variables));

  SLEQP_CALL(sleqp_vec_create_empty(&data->hessian_prod_func, num_variables));

  SLEQP_CALL(sleqp_vec_create_empty(&data->hessian_prod_cache, num_variables));

  SLEQP_CALL(sleqp_mat_create(&data->hessian_cons_prods,
                              num_variables,
                              num_constraints,
                              0));

  SLEQP_CALL(sleqp_vec_create_empty(&data->cons_grad_iterate, num_variables));

  SLEQP_CALL(
    sleqp_vec_create_empty(&data->cons_grad_check_iterate, num_variables));

  SLEQP_CALL(
    sleqp_vec_create_empty(&data->transposed_jacobian_product, num_variables));

  SLEQP_CALL(
    sleqp_vec_create_empty(&data->combined_cons_grad_iterate, num_variables));

  SLEQP_CALL(sleqp_vec_create_empty(&data->combined_cons_grad_check_iterate,
                                    num_variables));

  SLEQP_CALL(sleqp_vec_create(&data->multipliers, num_constraints, 1));

  SLEQP_CALL(sleqp_iterate_create(&data->check_iterate,
                                  problem,
                                  sleqp_problem_vars_lb(problem)));

  SLEQP_CALL(sleqp_vec_create_empty(&data->jac_row, num_variables));

  SLEQP_CALL(sleqp_vec_create_empty(&data->check_jac_row, num_variables));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_deriv_check_perform(SleqpDerivChecker* deriv_checker,
                          SleqpIterate* iterate,
                          SLEQP_DERIV_CHECK flags)
{
  if (flags == SLEQP_DERIV_CHECK_SKIP)
  {
    return SLEQP_OKAY;
  }

  deriv_checker->valid_deriv = true;

  SleqpProblem* problem = deriv_checker->problem;

  const int num_variables = sleqp_problem_num_vars(problem);

  {
    const bool first_order  = flags & SLEQP_DERIV_CHECK_FIRST;
    const bool second_order = flags
                              & (SLEQP_DERIV_CHECK_SECOND_SIMPLE
                                 | SLEQP_DERIV_CHECK_SECOND_EXHAUSTIVE);

    if (first_order && second_order)
    {
      sleqp_log_debug("Checking first and second order derivatives");
    }
    else
    {
      if (first_order)
      {
        sleqp_log_debug("Checking first order derivatives");
      }
      else if (second_order)
      {
        sleqp_log_debug("Checking second order derivatives");
      }
    }
  }

  for (int j = 0; j < num_variables; ++j)
  {
    SLEQP_CALL(eval_and_check_deriv(deriv_checker, iterate, flags, j));
  }

  if (!deriv_checker->valid_deriv)
  {
    sleqp_raise(SLEQP_INVALID_DERIV, "Invalid derivative");
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_deriv_checker_free(SleqpDerivChecker** star)
{
  SleqpDerivChecker* deriv_checker = *star;

  if (!deriv_checker)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_vec_free(&deriv_checker->check_jac_row));

  SLEQP_CALL(sleqp_vec_free(&deriv_checker->jac_row));

  SLEQP_CALL(sleqp_iterate_release(&deriv_checker->check_iterate));

  SLEQP_CALL(sleqp_vec_free(&deriv_checker->multipliers));

  SLEQP_CALL(sleqp_vec_free(&deriv_checker->combined_cons_grad_check_iterate));

  SLEQP_CALL(sleqp_vec_free(&deriv_checker->combined_cons_grad_iterate));

  SLEQP_CALL(sleqp_vec_free(&deriv_checker->transposed_jacobian_product));

  SLEQP_CALL(sleqp_vec_free(&deriv_checker->cons_grad_check_iterate));

  SLEQP_CALL(sleqp_vec_free(&deriv_checker->cons_grad_iterate));

  SLEQP_CALL(sleqp_mat_release(&deriv_checker->hessian_cons_prods));

  SLEQP_CALL(sleqp_vec_free(&deriv_checker->hessian_prod_cache));

  SLEQP_CALL(sleqp_vec_free(&deriv_checker->hessian_prod_func));

  SLEQP_CALL(sleqp_vec_free(&deriv_checker->hessian_prod_cons));

  SLEQP_CALL(sleqp_vec_free(&deriv_checker->hessian_prod));

  SLEQP_CALL(sleqp_vec_free(&deriv_checker->hessian_right));

  SLEQP_CALL(sleqp_vec_free(&deriv_checker->hessian_left));

  SLEQP_CALL(sleqp_vec_free(&deriv_checker->hessian_estimate));

  SLEQP_CALL(sleqp_vec_free(&deriv_checker->unit_direction));

  SLEQP_CALL(sleqp_settings_release(&deriv_checker->settings));

  SLEQP_CALL(sleqp_problem_release(&deriv_checker->problem));

  sleqp_free(star);

  return SLEQP_OKAY;
}
