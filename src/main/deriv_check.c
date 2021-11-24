#include "deriv_check.h"

#include "cmp.h"
#include "log.h"
#include "mem.h"
#include "problem.h"
#include "sparse/sparse_matrix.h"

struct SleqpDerivChecker
{
  SleqpProblem* problem;
  SleqpParams* params;

  SleqpSparseVec* unit_direction;

  SleqpSparseVec* hessian_estimate;

  SleqpSparseVec* hessian_left;
  SleqpSparseVec* hessian_right;

  SleqpSparseVec* hessian_prod;
  SleqpSparseVec* hessian_prod_cons;
  SleqpSparseVec* hessian_prod_func;

  SleqpSparseVec* hessian_prod_cache;

  SleqpSparseMatrix* hessian_cons_prods;

  SleqpSparseVec* cons_grad_iterate;
  SleqpSparseVec* cons_grad_check_iterate;

  SleqpSparseVec* transposed_jacobian_product;
  SleqpSparseVec* combined_cons_grad_iterate;
  SleqpSparseVec* combined_cons_grad_check_iterate;

  SleqpSparseVec* multipliers;

  bool valid_deriv;

  SleqpIterate* iterate;
  SleqpIterate* check_iterate;

  SleqpSparseVec* jac_row;
  SleqpSparseVec* check_jac_row;
};

static const double one = 1.;

static SLEQP_RETCODE
restore_iterate(SleqpDerivChecker* deriv_checker, SleqpIterate* iterate)
{
  SleqpProblem* problem = deriv_checker->problem;

  bool reject = false;

  int obl_grad_nnz = 0;
  int cons_val_nnz = 0;
  int cons_jac_nnz = 0;

  SLEQP_CALL(sleqp_problem_set_value(problem,
                                     sleqp_iterate_primal(iterate),
                                     SLEQP_VALUE_REASON_CHECKING_DERIV,
                                     &reject,
                                     &obl_grad_nnz,
                                     &cons_val_nnz,
                                     &cons_jac_nnz));

  if (reject)
  {
    sleqp_log_error("Function rejected restoration after derivative check");
    return SLEQP_INTERNAL_ERROR;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
set_check_iterate(SleqpDerivChecker* deriv_checker, bool* reject)
{
  SleqpIterate* check_iterate = deriv_checker->check_iterate;
  SleqpProblem* problem       = deriv_checker->problem;

  int obj_grad_nnz = 0;
  int cons_val_nnz = 0;
  int cons_jac_nnz = 0;

  SLEQP_CALL(sleqp_problem_set_value(problem,
                                     sleqp_iterate_primal(check_iterate),
                                     SLEQP_VALUE_REASON_CHECKING_DERIV,
                                     reject,
                                     &obj_grad_nnz,
                                     &cons_val_nnz,
                                     &cons_jac_nnz));

  SLEQP_CALL(sleqp_sparse_vector_reserve(sleqp_iterate_obj_grad(check_iterate),
                                         obj_grad_nnz));

  SLEQP_CALL(sleqp_sparse_vector_reserve(sleqp_iterate_cons_val(check_iterate),
                                         cons_val_nnz));

  SLEQP_CALL(sleqp_sparse_matrix_reserve(sleqp_iterate_cons_jac(check_iterate),
                                         cons_jac_nnz));

  return SLEQP_OKAY;
}

static double
get_perturbation(const SleqpDerivChecker* deriv_checker,
                 const SleqpIterate* iterate,
                 int j)
{
  SleqpParams* params = deriv_checker->params;

  SleqpSparseVec* primal = sleqp_iterate_primal(iterate);

  double base_perturbation
    = sleqp_params_value(params, SLEQP_PARAM_DERIV_PERTURBATION);

  double value = sleqp_sparse_vector_value_at(primal, j);

  value = SLEQP_ABS(value);

  value = SLEQP_MAX(value, 1.);

  return value * base_perturbation;
}

static SLEQP_RETCODE
create_unit_direction(SleqpSparseVec* direction, const int j)
{
  SLEQP_CALL(sleqp_sparse_vector_clear(direction));

  SLEQP_CALL(sleqp_sparse_vector_push(direction, j, 1.));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
create_perturbed_unit_direction(SleqpDerivChecker* deriv_checker,
                                const SleqpIterate* iterate,
                                const int j,
                                double* perturbation)
{
  SleqpSparseVec* unit_direction = deriv_checker->unit_direction;

  *perturbation = get_perturbation(deriv_checker, iterate, j);

  SLEQP_CALL(sleqp_sparse_vector_clear(unit_direction));

  SLEQP_CALL(sleqp_sparse_vector_push(unit_direction, j, *perturbation));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_hessian_cons_products(SleqpDerivChecker* deriv_checker, int j)
{
  SleqpProblem* problem = deriv_checker->problem;

  SleqpSparseVec* unit_direction = deriv_checker->unit_direction;

  const double zero_eps
    = sleqp_params_value(deriv_checker->params, SLEQP_PARAM_ZERO_EPS);

  SLEQP_CALL(create_unit_direction(unit_direction, j));

  SleqpSparseMatrix* hessian_prod_matrix = deriv_checker->hessian_cons_prods;

  SLEQP_CALL(sleqp_sparse_matrix_clear(hessian_prod_matrix));

  const int num_constraints = sleqp_problem_num_cons(problem);

  for (int k = 0; k < num_constraints; ++k)
  {
    SLEQP_CALL(create_unit_direction(deriv_checker->multipliers, k));

    SLEQP_CALL(sleqp_problem_hess_prod(problem,
                                       &one,
                                       unit_direction,
                                       deriv_checker->multipliers,
                                       deriv_checker->hessian_prod_cache));

    SLEQP_CALL(
      sleqp_sparse_vector_add_scaled(deriv_checker->hessian_prod_cache,
                                     deriv_checker->hessian_prod_func,
                                     1.,
                                     -1.,
                                     zero_eps,
                                     deriv_checker->hessian_prod_cons));

    const int nnz = sleqp_sparse_matrix_nnz(hessian_prod_matrix);

    SLEQP_CALL(
      sleqp_sparse_matrix_reserve(hessian_prod_matrix,
                                  nnz + deriv_checker->hessian_prod_cons->nnz));

    SLEQP_CALL(sleqp_sparse_matrix_push_column(hessian_prod_matrix, k));

    SLEQP_CALL(sleqp_sparse_matrix_push_vec(hessian_prod_matrix,
                                            k,
                                            deriv_checker->hessian_prod_cons));
  }

  assert(sleqp_sparse_matrix_is_valid(hessian_prod_matrix));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_hessian_products(SleqpDerivChecker* deriv_checker,
                         SLEQP_DERIV_CHECK flags,
                         SleqpIterate* iterate,
                         const int j)
{
  SleqpProblem* problem = deriv_checker->problem;

  SleqpSparseVec* unit_direction = deriv_checker->unit_direction;

  SLEQP_CALL(create_unit_direction(unit_direction, j));

  if (flags & SLEQP_DERIV_CHECK_SECOND_CONS
      || flags & SLEQP_DERIV_CHECK_SECOND_OBJ)
  {
    SLEQP_CALL(sleqp_sparse_vector_clear(deriv_checker->multipliers));

    SLEQP_CALL(sleqp_problem_hess_prod(problem,
                                       &one,
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
    SleqpSparseVec* multipliers = sleqp_iterate_cons_dual(iterate);

    SLEQP_CALL(sleqp_problem_hess_prod(problem,
                                       &one,
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
  SleqpIterate* check_iterate    = deriv_checker->check_iterate;
  SleqpSparseVec* unit_direction = deriv_checker->unit_direction;

  const double zero_eps
    = sleqp_params_value(deriv_checker->params, SLEQP_PARAM_ZERO_EPS);

  SLEQP_CALL(
    create_perturbed_unit_direction(deriv_checker, iterate, j, perturbation));

  SLEQP_CALL(sleqp_sparse_vector_add(sleqp_iterate_primal(iterate),
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

  double* func_ptr            = NULL;
  SleqpSparseVec* cons_val    = NULL;
  SleqpSparseVec* obj_grad    = NULL;
  SleqpSparseMatrix* cons_jac = NULL;

  if (flags & SLEQP_DERIV_CHECK_FIRST_OBJ)
  {
    func_ptr = &obj_val;
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
    sleqp_problem_eval(problem, NULL, func_ptr, obj_grad, cons_val, cons_jac));

  if (func_ptr)
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
    = sleqp_params_value(deriv_checker->params, SLEQP_PARAM_DERIV_TOL);

  const double func_diff
    = sleqp_iterate_obj_val(check_iterate) - sleqp_iterate_obj_val(iterate);

  const double actual_value = func_diff / perturbation;

  const double expected_value
    = sleqp_sparse_vector_value_at(sleqp_iterate_obj_grad(iterate), j);

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
  SleqpSparseMatrix* cons_jac = sleqp_iterate_cons_jac(iterate);

  const int num_constraints = sleqp_problem_num_cons(problem);

  const double tolerance
    = sleqp_params_value(deriv_checker->params, SLEQP_PARAM_DERIV_TOL);

  for (int i = 0; i < num_constraints; ++i)
  {
    const double expected_value = sleqp_sparse_matrix_value_at(cons_jac, i, j);

    const double lower_value
      = sleqp_sparse_vector_value_at(sleqp_iterate_cons_val(iterate), i);

    const double upper_value
      = sleqp_sparse_vector_value_at(sleqp_iterate_cons_val(check_iterate), i);

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
    = sleqp_params_value(deriv_checker->params, SLEQP_PARAM_DERIV_TOL);

  SLEQP_CALL(
    sleqp_sparse_vector_add_scaled(sleqp_iterate_obj_grad(iterate),
                                   sleqp_iterate_obj_grad(check_iterate),
                                   -1.,
                                   1.,
                                   0,
                                   deriv_checker->hessian_estimate));

  SLEQP_CALL(sleqp_sparse_vector_scale(deriv_checker->hessian_estimate,
                                       1. / perturbation));

  for (int k = 0; k < num_variables; ++k)
  {
    SLEQP_CALL(create_unit_direction(deriv_checker->hessian_left, k));

    double expected_value;

    SLEQP_CALL(sleqp_sparse_vector_dot(deriv_checker->hessian_left,
                                       deriv_checker->hessian_prod_func,
                                       &expected_value));

    const double actual_value
      = sleqp_sparse_vector_value_at(deriv_checker->hessian_estimate, k);

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
    = sleqp_params_value(deriv_checker->params, SLEQP_PARAM_ZERO_EPS);

  const double tolerance
    = sleqp_params_value(deriv_checker->params, SLEQP_PARAM_DERIV_TOL);

  SleqpSparseMatrix* cons_jac  = sleqp_iterate_cons_jac(iterate);
  SleqpSparseMatrix* check_jac = sleqp_iterate_cons_jac(check_iterate);

  SleqpSparseVec* expected_hessian_prod = deriv_checker->hessian_prod_cons;
  SleqpSparseVec* actual_hessian_prod   = deriv_checker->hessian_estimate;

  for (int i = 0; i < num_constraints; ++i)
  {
    SLEQP_CALL(sleqp_sparse_matrix_col(deriv_checker->hessian_cons_prods,
                                       i,
                                       expected_hessian_prod));

    SLEQP_CALL(create_unit_direction(deriv_checker->multipliers, i));

    SLEQP_CALL(
      sleqp_sparse_matrix_trans_vector_product(cons_jac,
                                               deriv_checker->multipliers,
                                               zero_eps,
                                               deriv_checker->jac_row));

    SLEQP_CALL(
      sleqp_sparse_matrix_trans_vector_product(check_jac,
                                               deriv_checker->multipliers,
                                               zero_eps,
                                               deriv_checker->check_jac_row));

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(deriv_checker->check_jac_row,
                                              deriv_checker->jac_row,
                                              1.,
                                              -1.,
                                              zero_eps,
                                              actual_hessian_prod));

    SLEQP_CALL(
      sleqp_sparse_vector_scale(actual_hessian_prod, 1. / perturbation));

    // validation
    for (int k = 0; k < num_variables; ++k)
    {
      const double expected_value
        = sleqp_sparse_vector_value_at(expected_hessian_prod, k);

      const double actual_value
        = sleqp_sparse_vector_value_at(actual_hessian_prod, k);

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
                           SleqpSparseVec* multipliers,
                           SleqpIterate* iterate,
                           SleqpSparseVec* result)
{
  SleqpParams* params = deriv_checker->params;

  const double zero_eps = sleqp_params_value(params, SLEQP_PARAM_ZERO_EPS);

  SLEQP_CALL(sleqp_sparse_matrix_trans_vector_product(
    sleqp_iterate_cons_jac(iterate),
    multipliers,
    zero_eps,
    deriv_checker->transposed_jacobian_product));

  SLEQP_CALL(sleqp_sparse_vector_add(sleqp_iterate_obj_grad(iterate),
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
  SleqpSparseVec* multipliers = sleqp_iterate_cons_dual(iterate);

  const int num_variables = sleqp_problem_num_vars(problem);

  const double tolerance
    = sleqp_params_value(deriv_checker->params, SLEQP_PARAM_DERIV_TOL);

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

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(
      deriv_checker->combined_cons_grad_iterate,
      deriv_checker->combined_cons_grad_check_iterate,
      -1.,
      1.,
      0,
      deriv_checker->hessian_estimate));

    SLEQP_CALL(sleqp_sparse_vector_scale(deriv_checker->hessian_estimate,
                                         1. / perturbation));
  }

  // Validation
  for (int k = 0; k < num_variables; ++k)
  {
    const double actual_value
      = sleqp_sparse_vector_value_at(deriv_checker->hessian_estimate, k);

    const double expected_value
      = sleqp_sparse_vector_value_at(deriv_checker->hessian_prod, k);

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
                           SleqpParams* params)
{
  SLEQP_CALL(sleqp_malloc(deriv_checker));

  SleqpDerivChecker* data = *deriv_checker;

  const int num_constraints = sleqp_problem_num_cons(problem);
  const int num_variables   = sleqp_problem_num_vars(problem);

  data->problem = problem;
  SLEQP_CALL(sleqp_problem_capture(data->problem));

  SLEQP_CALL(sleqp_params_capture(params));
  data->params = params;

  data->iterate = NULL;

  // this will be a unit vector
  SLEQP_CALL(
    sleqp_sparse_vector_create(&data->unit_direction, num_variables, 1));

  SLEQP_CALL(
    sleqp_sparse_vector_create_empty(&data->hessian_estimate, num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->hessian_left, num_variables, 1));

  SLEQP_CALL(
    sleqp_sparse_vector_create(&data->hessian_right, num_variables, 1));

  SLEQP_CALL(
    sleqp_sparse_vector_create_empty(&data->hessian_prod, num_variables));

  SLEQP_CALL(
    sleqp_sparse_vector_create_empty(&data->hessian_prod_cons, num_variables));

  SLEQP_CALL(
    sleqp_sparse_vector_create_empty(&data->hessian_prod_func, num_variables));

  SLEQP_CALL(
    sleqp_sparse_vector_create_empty(&data->hessian_prod_cache, num_variables));

  SLEQP_CALL(sleqp_sparse_matrix_create(&data->hessian_cons_prods,
                                        num_variables,
                                        num_constraints,
                                        0));

  SLEQP_CALL(
    sleqp_sparse_vector_create_empty(&data->cons_grad_iterate, num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->cons_grad_check_iterate,
                                              num_variables));

  SLEQP_CALL(
    sleqp_sparse_vector_create_empty(&data->transposed_jacobian_product,
                                     num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->combined_cons_grad_iterate,
                                              num_variables));

  SLEQP_CALL(
    sleqp_sparse_vector_create_empty(&data->combined_cons_grad_check_iterate,
                                     num_variables));

  SLEQP_CALL(
    sleqp_sparse_vector_create(&data->multipliers, num_constraints, 1));

  SLEQP_CALL(sleqp_iterate_create(&data->check_iterate,
                                  problem,
                                  sleqp_problem_vars_lb(problem)));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->jac_row, num_variables));

  SLEQP_CALL(
    sleqp_sparse_vector_create_empty(&data->check_jac_row, num_variables));

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

  for (int j = 0; j < num_variables; ++j)
  {
    SLEQP_CALL(eval_and_check_deriv(deriv_checker, iterate, flags, j));
  }

  if (!deriv_checker->valid_deriv)
  {
    return SLEQP_INVALID_DERIV;
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

  SLEQP_CALL(sleqp_sparse_vector_free(&deriv_checker->check_jac_row));

  SLEQP_CALL(sleqp_sparse_vector_free(&deriv_checker->jac_row));

  SLEQP_CALL(sleqp_iterate_release(&deriv_checker->check_iterate));

  SLEQP_CALL(sleqp_sparse_vector_free(&deriv_checker->multipliers));

  SLEQP_CALL(
    sleqp_sparse_vector_free(&deriv_checker->combined_cons_grad_check_iterate));

  SLEQP_CALL(
    sleqp_sparse_vector_free(&deriv_checker->combined_cons_grad_iterate));

  SLEQP_CALL(
    sleqp_sparse_vector_free(&deriv_checker->transposed_jacobian_product));

  SLEQP_CALL(sleqp_sparse_vector_free(&deriv_checker->cons_grad_check_iterate));

  SLEQP_CALL(sleqp_sparse_vector_free(&deriv_checker->cons_grad_iterate));

  SLEQP_CALL(sleqp_sparse_matrix_release(&deriv_checker->hessian_cons_prods));

  SLEQP_CALL(sleqp_sparse_vector_free(&deriv_checker->hessian_prod_cache));

  SLEQP_CALL(sleqp_sparse_vector_free(&deriv_checker->hessian_prod_func));

  SLEQP_CALL(sleqp_sparse_vector_free(&deriv_checker->hessian_prod_cons));

  SLEQP_CALL(sleqp_sparse_vector_free(&deriv_checker->hessian_prod));

  SLEQP_CALL(sleqp_sparse_vector_free(&deriv_checker->hessian_right));

  SLEQP_CALL(sleqp_sparse_vector_free(&deriv_checker->hessian_left));

  SLEQP_CALL(sleqp_sparse_vector_free(&deriv_checker->hessian_estimate));

  SLEQP_CALL(sleqp_sparse_vector_free(&deriv_checker->unit_direction));

  SLEQP_CALL(sleqp_params_release(&deriv_checker->params));

  SLEQP_CALL(sleqp_problem_release(&deriv_checker->problem));

  sleqp_free(star);

  return SLEQP_OKAY;
}
