#include "problem_scaling.h"

#include "cmp.h"
#include "dyn.h"
#include "error.h"
#include "func.h"
#include "log.h"
#include "lsq.h"
#include "math_error.h"
#include "mem.h"
#include "sparse/mat.h"

struct SleqpProblemScaling
{
  int refcount;

  SleqpScaling* scaling;
  SleqpProblem* problem;
  SleqpSettings* settings;
  SleqpFunc* func;

  SleqpFunc* scaled_func;
  SleqpProblem* scaled_problem;

  SleqpVec* unscaled_value;

  SleqpVec* scaled_direction;

  double* scaled_cons_weights;
  SleqpVec* scaled_cons_duals;
};

static SLEQP_RETCODE
scaled_func_set_value(SleqpFunc* func,
                      SleqpVec* scaled_value,
                      SLEQP_VALUE_REASON reason,
                      bool* reject,
                      void* func_data)
{
  SleqpProblemScaling* problem_scaling = (SleqpProblemScaling*)func_data;
  SleqpScaling* scaling                = problem_scaling->scaling;

  {
    const int error_flags
      = sleqp_settings_enum_value(problem_scaling->settings,
                                 SLEQP_SETTINGS_ENUM_FLOAT_ERROR_FLAGS);

    const int warn_flags
      = sleqp_settings_enum_value(problem_scaling->settings,
                                 SLEQP_SETTINGS_ENUM_FLOAT_WARNING_FLAGS);

    SLEQP_INIT_MATH_CHECK;

    SLEQP_CALL(sleqp_vec_copy(scaled_value, problem_scaling->unscaled_value));

    SLEQP_CALL(sleqp_unscale_point(scaling, problem_scaling->unscaled_value));

    SLEQP_MATH_CHECK(error_flags, warn_flags);
  }

  SLEQP_CALL(sleqp_func_set_value(problem_scaling->func,
                                  problem_scaling->unscaled_value,
                                  reason,
                                  reject));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
scaled_func_nonzeros(SleqpFunc* func,
                     int* obj_grad_nnz,
                     int* cons_val_nnz,
                     int* cons_jac_nnz,
                     int* hess_prod_nnz,
                     void* func_data)
{
  SleqpProblemScaling* problem_scaling = (SleqpProblemScaling*)func_data;

  SLEQP_CALL(sleqp_func_nonzeros(problem_scaling->func,
                                 obj_grad_nnz,
                                 cons_val_nnz,
                                 cons_jac_nnz,
                                 hess_prod_nnz));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
scaled_func_obj_val(SleqpFunc* func, double* obj_val, void* func_data)
{
  SleqpProblemScaling* problem_scaling = (SleqpProblemScaling*)func_data;
  SleqpScaling* scaling                = problem_scaling->scaling;

  SLEQP_CALL(sleqp_func_obj_val(problem_scaling->func, obj_val));

  (*obj_val) = sleqp_scale_obj_val(scaling, (*obj_val));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
scaled_func_obj_grad(SleqpFunc* func, SleqpVec* obj_grad, void* func_data)
{
  SleqpProblemScaling* problem_scaling = (SleqpProblemScaling*)func_data;
  SleqpScaling* scaling                = problem_scaling->scaling;

  SLEQP_CALL(sleqp_func_obj_grad(problem_scaling->func, obj_grad));

  SLEQP_CALL(sleqp_scale_obj_grad(scaling, obj_grad));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
scaled_func_cons_val(SleqpFunc* func, SleqpVec* cons_val, void* func_data)
{
  SleqpProblemScaling* problem_scaling = (SleqpProblemScaling*)func_data;
  SleqpScaling* scaling                = problem_scaling->scaling;

  SLEQP_CALL(sleqp_func_cons_val(problem_scaling->func, cons_val));

  SLEQP_CALL(sleqp_scale_cons_val(scaling, cons_val));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
scaled_func_cons_jac(SleqpFunc* func, SleqpMat* cons_jac, void* func_data)
{
  SleqpProblemScaling* problem_scaling = (SleqpProblemScaling*)func_data;
  SleqpScaling* scaling                = problem_scaling->scaling;

  SLEQP_CALL(sleqp_func_cons_jac(problem_scaling->func, cons_jac));

  SLEQP_CALL(sleqp_scale_cons_jac(scaling, cons_jac));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
scaled_func_hess_prod(SleqpFunc* func,
                      const SleqpVec* direction,
                      const SleqpVec* cons_duals,
                      SleqpVec* product,
                      void* func_data)
{
  SleqpProblemScaling* problem_scaling = (SleqpProblemScaling*)func_data;
  SleqpScaling* scaling                = problem_scaling->scaling;

  const int error_flags
    = sleqp_settings_enum_value(problem_scaling->settings,
                               SLEQP_SETTINGS_ENUM_FLOAT_ERROR_FLAGS);

  const int warn_flags
    = sleqp_settings_enum_value(problem_scaling->settings,
                               SLEQP_SETTINGS_ENUM_FLOAT_WARNING_FLAGS);

  SLEQP_CALL(sleqp_vec_copy(direction, problem_scaling->scaled_direction));

  SLEQP_CALL(sleqp_vec_copy(cons_duals, problem_scaling->scaled_cons_duals));

  {
    SLEQP_INIT_MATH_CHECK;

    SLEQP_CALL(
      sleqp_unscale_hessian_direction(scaling,
                                      problem_scaling->scaled_direction,
                                      problem_scaling->scaled_cons_duals));

    SLEQP_MATH_CHECK(error_flags, warn_flags);
  }

  SLEQP_CALL(sleqp_func_hess_prod(problem_scaling->func,
                                  problem_scaling->scaled_direction,
                                  problem_scaling->scaled_cons_duals,
                                  product));

  {
    SLEQP_INIT_MATH_CHECK;

    SLEQP_CALL(sleqp_scale_hessian_product(scaling, product));

    SLEQP_MATH_CHECK(error_flags, warn_flags);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
scaled_lsq_func_nonzeros(SleqpFunc* func,
                         int* residual_nnz,
                         int* jac_fwd_nnz,
                         int* jac_adj_nnz,
                         int* cons_val_nnz,
                         int* cons_jac_nnz,
                         void* func_data)
{
  SleqpProblemScaling* problem_scaling = (SleqpProblemScaling*)func_data;

  SLEQP_CALL(sleqp_lsq_func_nonzeros(problem_scaling->func,
                                     residual_nnz,
                                     jac_fwd_nnz,
                                     jac_adj_nnz,
                                     cons_val_nnz,
                                     cons_jac_nnz));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
scaled_lsq_func_residuals(SleqpFunc* func, SleqpVec* residuals, void* func_data)
{
  SleqpProblemScaling* problem_scaling = (SleqpProblemScaling*)func_data;
  SleqpScaling* scaling                = problem_scaling->scaling;

  SLEQP_CALL(sleqp_lsq_func_residuals(problem_scaling->func, residuals));

  SLEQP_CALL(sleqp_scale_lsq_residuals(scaling, residuals));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
scaled_lsq_func_jac_forward(SleqpFunc* func,
                            const SleqpVec* forward_direction,
                            SleqpVec* product,
                            void* func_data)
{
  SleqpProblemScaling* problem_scaling = (SleqpProblemScaling*)func_data;
  SleqpScaling* scaling                = problem_scaling->scaling;

  SleqpVec* scaled_direction = problem_scaling->scaled_direction;

  SLEQP_CALL(sleqp_vec_resize(scaled_direction, forward_direction->dim));

  SLEQP_CALL(sleqp_vec_copy(forward_direction, scaled_direction));

  SLEQP_CALL(sleqp_scale_lsq_forward_direction(scaling, scaled_direction));

  SLEQP_CALL(sleqp_lsq_func_jac_forward(problem_scaling->func,
                                        scaled_direction,
                                        product));

  SLEQP_CALL(sleqp_scale_lsq_adjoint_direction(scaling, product));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
scaled_lsq_func_jac_adjoint(SleqpFunc* func,
                            const SleqpVec* adjoint_direction,
                            SleqpVec* product,
                            void* func_data)
{
  SleqpProblemScaling* problem_scaling = (SleqpProblemScaling*)func_data;
  SleqpScaling* scaling                = problem_scaling->scaling;

  SleqpVec* scaled_direction = problem_scaling->scaled_direction;

  SLEQP_CALL(sleqp_vec_resize(scaled_direction, adjoint_direction->dim));

  SLEQP_CALL(sleqp_vec_copy(adjoint_direction, scaled_direction));

  SLEQP_CALL(sleqp_scale_lsq_adjoint_direction(scaling, scaled_direction));

  SLEQP_CALL(sleqp_lsq_func_jac_adjoint(problem_scaling->func,
                                        scaled_direction,
                                        product));

  SLEQP_CALL(sleqp_scale_lsq_forward_direction(scaling, product));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
scaled_dyn_func_eval(SleqpFunc* func,
                     double* obj_val,
                     SleqpVec* cons_val,
                     double* error,
                     void* func_data)
{
  SleqpProblemScaling* problem_scaling = (SleqpProblemScaling*)func_data;
  SleqpScaling* scaling                = problem_scaling->scaling;

  SLEQP_CALL(
    sleqp_dyn_func_eval(problem_scaling->func, obj_val, cons_val, error));

  SLEQP_CALL(sleqp_scale_cons_val(scaling, cons_val));
  (*obj_val) = sleqp_scale_obj_val(scaling, (*obj_val));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
scaled_dyn_func_set_error_bound(SleqpFunc* func,
                                double error_bound,
                                void* func_data)
{
  SleqpProblemScaling* problem_scaling = (SleqpProblemScaling*)func_data;

  SLEQP_CALL(
    sleqp_dyn_func_set_error_bound(problem_scaling->func, error_bound));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
scaled_dyn_func_set_obj_weight(SleqpFunc* func,
                               double obj_weight,
                               void* func_data)
{
  SleqpProblemScaling* problem_scaling = (SleqpProblemScaling*)func_data;
  SleqpScaling* scaling                = problem_scaling->scaling;

  const double scaled_obj_weight = sleqp_scale_obj_weight(scaling, obj_weight);

  SLEQP_CALL(
    sleqp_dyn_func_set_obj_weight(problem_scaling->func, scaled_obj_weight));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
scaled_dyn_func_set_cons_weights(SleqpFunc* func,
                                 const double* cons_weights,
                                 void* func_data)
{
  SleqpProblemScaling* problem_scaling = (SleqpProblemScaling*)func_data;
  SleqpScaling* scaling                = problem_scaling->scaling;
  double* scaled_cons_weights          = problem_scaling->scaled_cons_weights;

  SLEQP_CALL(
    sleqp_scale_cons_weights(scaling, cons_weights, scaled_cons_weights));

  SLEQP_CALL(sleqp_dyn_func_set_cons_weights(problem_scaling->func,
                                             scaled_cons_weights));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
func_create(SleqpProblemScaling* problem_scaling)
{
  SleqpProblem* problem = problem_scaling->problem;

  const int num_variables   = sleqp_problem_num_vars(problem);
  const int num_constraints = sleqp_problem_num_cons(problem);

  SleqpFuncCallbacks callbacks = {.set_value = scaled_func_set_value,
                                  .nonzeros  = scaled_func_nonzeros,
                                  .obj_val   = scaled_func_obj_val,
                                  .obj_grad  = scaled_func_obj_grad,
                                  .cons_val  = scaled_func_cons_val,
                                  .cons_jac  = scaled_func_cons_jac,
                                  .hess_prod = scaled_func_hess_prod,
                                  .func_free = NULL};

  SLEQP_CALL(sleqp_func_create(&(problem_scaling->scaled_func),
                               &callbacks,
                               num_variables,
                               num_constraints,
                               problem_scaling));

  SLEQP_CALL(sleqp_hess_struct_copy(
    sleqp_func_hess_struct(problem_scaling->func),
    sleqp_func_hess_struct(problem_scaling->scaled_func)));

  SLEQP_CALL(
    sleqp_func_flags_copy(problem_scaling->func,
                          problem_scaling->scaled_func,
                          SLEQP_FUNC_HESS_INEXACT | SLEQP_FUNC_HESS_PSD));

  SLEQP_CALL(
    sleqp_func_flags_add(problem_scaling->scaled_func, SLEQP_FUNC_INTERNAL));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
lsq_func_create(SleqpProblemScaling* problem_scaling)
{
  SleqpProblem* problem = problem_scaling->problem;

  const int num_variables   = sleqp_problem_num_vars(problem);
  const int num_constraints = sleqp_problem_num_cons(problem);

  SleqpLSQCallbacks callbacks = {.set_value       = scaled_func_set_value,
                                 .lsq_nonzeros    = scaled_lsq_func_nonzeros,
                                 .lsq_residuals   = scaled_lsq_func_residuals,
                                 .lsq_jac_forward = scaled_lsq_func_jac_forward,
                                 .lsq_jac_adjoint = scaled_lsq_func_jac_adjoint,
                                 .cons_val        = scaled_func_cons_val,
                                 .cons_jac        = scaled_func_cons_jac,
                                 .func_free       = NULL};

  const double levenberg_marquardt
    = sleqp_lsq_func_get_levenberg_marquardt(problem_scaling->func);
  const int num_residuals = sleqp_lsq_func_num_residuals(problem_scaling->func);

  SLEQP_CALL(sleqp_lsq_func_create(&(problem_scaling->scaled_func),
                                   &callbacks,
                                   num_variables,
                                   num_constraints,
                                   num_residuals,
                                   levenberg_marquardt,
                                   problem_scaling->settings,
                                   problem_scaling));

  SLEQP_CALL(
    sleqp_func_flags_add(problem_scaling->scaled_func, SLEQP_FUNC_INTERNAL));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
dyn_func_create(SleqpProblemScaling* problem_scaling)
{
  SleqpProblem* problem = problem_scaling->problem;

  const int num_variables   = sleqp_problem_num_vars(problem);
  const int num_constraints = sleqp_problem_num_cons(problem);

  SleqpDynFuncCallbacks callbacks
    = {.set_value        = scaled_func_set_value,
       .nonzeros         = scaled_func_nonzeros,
       .set_error_bound  = scaled_dyn_func_set_error_bound,
       .set_obj_weight   = scaled_dyn_func_set_obj_weight,
       .set_cons_weights = scaled_dyn_func_set_cons_weights,
       .eval             = scaled_dyn_func_eval,
       .obj_grad         = scaled_func_obj_grad,
       .cons_jac         = scaled_func_cons_jac,
       .hess_prod        = scaled_func_hess_prod,
       .func_free        = NULL};

  SLEQP_CALL(sleqp_dyn_func_create(&(problem_scaling->scaled_func),
                                   &callbacks,
                                   num_variables,
                                   num_constraints,
                                   problem_scaling));

  SLEQP_CALL(
    sleqp_func_flags_add(problem_scaling->scaled_func, SLEQP_FUNC_INTERNAL));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_problem_scaling_create(SleqpProblemScaling** star,
                             SleqpScaling* scaling,
                             SleqpProblem* problem,
                             SleqpSettings* settings)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpProblemScaling* problem_scaling = *star;

  *problem_scaling = (SleqpProblemScaling){0};

  const int num_variables   = sleqp_problem_num_vars(problem);
  const int num_constraints = sleqp_problem_num_cons(problem);

  if (num_variables != sleqp_scaling_num_vars(scaling))
  {
    sleqp_raise(SLEQP_ILLEGAL_ARGUMENT,
                "Invalid number of variables provided to scaled problem, "
                "expected %d, actual: %d",
                num_variables,
                sleqp_scaling_num_vars(scaling));
  }

  if (num_constraints != sleqp_scaling_num_cons(scaling))
  {
    sleqp_raise(SLEQP_ILLEGAL_ARGUMENT,
                "Invalid number of constraints provided to scaled problem, "
                "expected %d, actual %d",
                num_constraints,
                sleqp_scaling_num_cons(scaling));
  }

  problem_scaling->refcount = 1;

  problem_scaling->problem = problem;
  SLEQP_CALL(sleqp_problem_capture(problem_scaling->problem));

  problem_scaling->func = sleqp_problem_func(problem);

  SLEQP_CALL(sleqp_settings_capture(settings));
  problem_scaling->settings = settings;

  problem_scaling->scaling = scaling;

  SLEQP_CALL(sleqp_scaling_capture(problem_scaling->scaling));

  switch (sleqp_func_get_type(problem_scaling->func))
  {
  case SLEQP_FUNC_TYPE_REGULAR:
    SLEQP_CALL(func_create(problem_scaling));
    break;
  case SLEQP_FUNC_TYPE_LSQ:
    SLEQP_CALL(lsq_func_create(problem_scaling));
    break;
  case SLEQP_FUNC_TYPE_DYNAMIC:
    SLEQP_CALL(dyn_func_create(problem_scaling));
    SLEQP_CALL(sleqp_alloc_array(&problem_scaling->scaled_cons_weights,
                                 num_constraints));
    break;
  }

  SLEQP_CALL(sleqp_problem_create(&(problem_scaling->scaled_problem),
                                  problem_scaling->scaled_func,
                                  sleqp_problem_vars_lb(problem),
                                  sleqp_problem_vars_ub(problem),
                                  sleqp_problem_general_lb(problem),
                                  sleqp_problem_general_ub(problem),
                                  sleqp_problem_linear_coeffs(problem),
                                  sleqp_problem_linear_lb(problem),
                                  sleqp_problem_linear_ub(problem),
                                  problem_scaling->settings));

  SLEQP_CALL(
    sleqp_vec_create_empty(&(problem_scaling->unscaled_value), num_variables));

  SLEQP_CALL(sleqp_vec_create_empty(&(problem_scaling->scaled_direction),
                                    num_variables));

  SLEQP_CALL(sleqp_vec_create_empty(&(problem_scaling->scaled_cons_duals),
                                    num_constraints));

  return SLEQP_OKAY;
}

SleqpProblem*
sleqp_problem_scaling_get_problem(SleqpProblemScaling* problem_scaling)
{
  return problem_scaling->scaled_problem;
}

SLEQP_RETCODE
sleqp_problem_scaling_flush(SleqpProblemScaling* problem_scaling)
{
  SleqpProblem* problem        = problem_scaling->problem;
  SleqpScaling* scaling        = problem_scaling->scaling;
  SleqpProblem* scaled_problem = problem_scaling->scaled_problem;

  const int error_flags
    = sleqp_settings_enum_value(problem_scaling->settings,
                               SLEQP_SETTINGS_ENUM_FLOAT_ERROR_FLAGS);

  const int warn_flags
    = sleqp_settings_enum_value(problem_scaling->settings,
                               SLEQP_SETTINGS_ENUM_FLOAT_WARNING_FLAGS);

  SLEQP_INIT_MATH_CHECK;

  SLEQP_CALL(sleqp_vec_copy(sleqp_problem_vars_lb(problem),
                            sleqp_problem_vars_lb(scaled_problem)));

  SLEQP_CALL(sleqp_scale_point(scaling, sleqp_problem_vars_lb(scaled_problem)));

  SLEQP_CALL(sleqp_vec_copy(sleqp_problem_vars_ub(problem),
                            sleqp_problem_vars_ub(scaled_problem)));

  SLEQP_CALL(sleqp_scale_point(scaling, sleqp_problem_vars_ub(scaled_problem)));

  SLEQP_CALL(sleqp_vec_copy(sleqp_problem_cons_lb(problem),
                            sleqp_problem_cons_lb(scaled_problem)));

  SLEQP_CALL(
    sleqp_scale_cons_val(scaling, sleqp_problem_cons_lb(scaled_problem)));

  SLEQP_CALL(sleqp_vec_copy(sleqp_problem_cons_ub(problem),
                            sleqp_problem_cons_ub(scaled_problem)));

  SLEQP_CALL(
    sleqp_scale_cons_val(scaling, sleqp_problem_cons_ub(scaled_problem)));

  SLEQP_CALL(sleqp_vec_copy(sleqp_problem_general_lb(problem),
                            sleqp_problem_general_lb(scaled_problem)));

  SLEQP_CALL(
    sleqp_scale_cons_general(scaling,
                             sleqp_problem_general_lb(scaled_problem)));

  SLEQP_CALL(sleqp_vec_copy(sleqp_problem_general_ub(problem),
                            sleqp_problem_general_ub(scaled_problem)));

  SLEQP_CALL(
    sleqp_scale_cons_general(scaling,
                             sleqp_problem_general_ub(scaled_problem)));

  SLEQP_CALL(sleqp_vec_copy(sleqp_problem_linear_lb(problem),
                            sleqp_problem_linear_lb(scaled_problem)));

  SLEQP_CALL(
    sleqp_scale_cons_linear(scaling, sleqp_problem_linear_lb(scaled_problem)));

  SLEQP_CALL(sleqp_vec_copy(sleqp_problem_linear_ub(problem),
                            sleqp_problem_linear_ub(scaled_problem)));

  SLEQP_CALL(
    sleqp_scale_cons_linear(scaling, sleqp_problem_linear_ub(scaled_problem)));

  SLEQP_CALL(sleqp_mat_copy(sleqp_problem_linear_coeffs(problem),
                            sleqp_problem_linear_coeffs(scaled_problem)));

  SLEQP_CALL(
    sleqp_scale_linear_coeffs(scaling,
                              sleqp_problem_linear_coeffs(scaled_problem)));

  SLEQP_MATH_CHECK(error_flags, warn_flags);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
problem_scaling_free(SleqpProblemScaling** star)
{
  SleqpProblemScaling* problem_scaling = *star;

  if (!problem_scaling)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_vec_free(&(problem_scaling->scaled_cons_duals)));

  sleqp_free(&problem_scaling->scaled_cons_weights);

  SLEQP_CALL(sleqp_vec_free(&(problem_scaling->scaled_direction)));

  SLEQP_CALL(sleqp_vec_free(&(problem_scaling->unscaled_value)));

  SLEQP_CALL(sleqp_problem_release(&(problem_scaling->scaled_problem)));

  SLEQP_CALL(sleqp_func_release(&(problem_scaling->scaled_func)));

  SLEQP_CALL(sleqp_scaling_release(&problem_scaling->scaling));

  SLEQP_CALL(sleqp_settings_release(&problem_scaling->settings));

  SLEQP_CALL(sleqp_problem_release(&problem_scaling->problem));

  sleqp_free(star);

  *star = NULL;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_problem_scaling_capture(SleqpProblemScaling* scaling)
{
  ++scaling->refcount;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_problem_scaling_release(SleqpProblemScaling** star)
{
  SleqpProblemScaling* problem_scaling = *star;

  if (!problem_scaling)
  {
    return SLEQP_OKAY;
  }

  if (--problem_scaling->refcount == 0)
  {
    SLEQP_CALL(problem_scaling_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
