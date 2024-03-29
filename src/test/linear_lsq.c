#include "linear_lsq.h"

#include "cmp.h"
#include "lsq.h"
#include "sparse/mat.h"

#define LINEAR_LSQ_NUM_VARIABLES 2
#define LINEAR_LSQ_NUM_RESIDUALS 3

const int linear_lsq_num_variables   = LINEAR_LSQ_NUM_VARIABLES;
const int linear_lsq_num_constraints = 0;
const int linear_lsq_num_residuals   = LINEAR_LSQ_NUM_RESIDUALS;

SleqpMat* linear_lsq_matrix;
SleqpVec* linear_lsq_rhs;

SleqpSettings* linear_lsq_settings;
SleqpFunc* linear_lsq_func;

SleqpVec* linear_lsq_var_lb;
SleqpVec* linear_lsq_var_ub;
SleqpVec* linear_lsq_cons_lb;
SleqpVec* linear_lsq_cons_ub;
SleqpVec* linear_lsq_initial;
SleqpVec* linear_lsq_optimal;

static SleqpVec* linear_lsq_current;
static SleqpVec* linear_lsq_forward;

double dense_cache_forward[LINEAR_LSQ_NUM_RESIDUALS];

static SLEQP_RETCODE
func_set(SleqpFunc* func,
         SleqpVec* value,
         SLEQP_VALUE_REASON reason,
         bool* reject,
         void* func_data)
{
  SLEQP_CALL(sleqp_vec_copy(value, linear_lsq_current));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
lsq_residuals(SleqpFunc* func, SleqpVec* residual, void* func_data)
{
  const double zero_eps
    = sleqp_settings_real_value(linear_lsq_settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  SLEQP_CALL(sleqp_mat_mult_vec(linear_lsq_matrix,
                                linear_lsq_current,
                                dense_cache_forward));

  SLEQP_CALL(sleqp_vec_set_from_raw(linear_lsq_forward,
                                    dense_cache_forward,
                                    linear_lsq_num_residuals,
                                    zero_eps));

  SLEQP_CALL(sleqp_vec_add_scaled(linear_lsq_forward,
                                  linear_lsq_rhs,
                                  1.,
                                  -1.,
                                  zero_eps,
                                  residual));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
lsq_jac_forward(SleqpFunc* func,
                const SleqpVec* forward_direction,
                SleqpVec* product,
                void* func_data)
{
  const double zero_eps
    = sleqp_settings_real_value(linear_lsq_settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  SLEQP_CALL(sleqp_mat_mult_vec(linear_lsq_matrix,
                                forward_direction,
                                dense_cache_forward));

  SLEQP_CALL(sleqp_vec_set_from_raw(product,
                                    dense_cache_forward,
                                    linear_lsq_num_residuals,
                                    zero_eps));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
lsq_jac_adjoint(SleqpFunc* func,
                const SleqpVec* adjoint_direction,
                SleqpVec* product,
                void* func_data)
{
  const double zero_eps
    = sleqp_settings_real_value(linear_lsq_settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  SLEQP_CALL(sleqp_mat_mult_vec_trans(linear_lsq_matrix,
                                      adjoint_direction,
                                      zero_eps,
                                      product));

  return SLEQP_OKAY;
}

void
linear_lsq_setup()
{
  const double inf = sleqp_infinity();

  ASSERT_CALL(
    sleqp_mat_create(&linear_lsq_matrix,
                     linear_lsq_num_residuals,
                     linear_lsq_num_variables,
                     linear_lsq_num_residuals * linear_lsq_num_variables));

  {
    ASSERT_CALL(sleqp_mat_clear(linear_lsq_matrix));

    ASSERT_CALL(sleqp_mat_push(linear_lsq_matrix, 0, 0, 1.));
    ASSERT_CALL(sleqp_mat_push(linear_lsq_matrix, 1, 0, 1.));
    ASSERT_CALL(sleqp_mat_push(linear_lsq_matrix, 2, 0, 1.));

    ASSERT_CALL(sleqp_mat_push_col(linear_lsq_matrix, 1));

    ASSERT_CALL(sleqp_mat_push(linear_lsq_matrix, 1, 1, 1.));
    ASSERT_CALL(sleqp_mat_push(linear_lsq_matrix, 2, 1, 2.));
  }

  ASSERT_CALL(sleqp_vec_create_full(&linear_lsq_rhs, linear_lsq_num_residuals));

  {
    double values[] = {6., 0., 0.};

    ASSERT_CALL(sleqp_vec_set_from_raw(linear_lsq_rhs,
                                       values,
                                       linear_lsq_num_residuals,
                                       0.));
  }

  ASSERT_CALL(
    sleqp_vec_create_full(&linear_lsq_var_lb, linear_lsq_num_variables));

  {
    double values[] = {-inf, -inf, -inf};

    ASSERT_CALL(sleqp_vec_set_from_raw(linear_lsq_var_lb,
                                       values,
                                       linear_lsq_num_variables,
                                       0.));
  }

  ASSERT_CALL(
    sleqp_vec_create_full(&linear_lsq_var_ub, linear_lsq_num_variables));

  {
    double values[] = {inf, inf, inf};

    ASSERT_CALL(sleqp_vec_set_from_raw(linear_lsq_var_ub,
                                       values,
                                       linear_lsq_num_variables,
                                       0.));
  }

  ASSERT_CALL(
    sleqp_vec_create_empty(&linear_lsq_cons_lb, linear_lsq_num_constraints));

  ASSERT_CALL(
    sleqp_vec_create_empty(&linear_lsq_cons_ub, linear_lsq_num_constraints));

  ASSERT_CALL(
    sleqp_vec_create_empty(&linear_lsq_initial, linear_lsq_num_variables));

  ASSERT_CALL(
    sleqp_vec_create_empty(&linear_lsq_optimal, linear_lsq_num_variables));

  {
    double values[] = {5., -3.};

    ASSERT_CALL(sleqp_vec_set_from_raw(linear_lsq_optimal,
                                       values,
                                       linear_lsq_num_variables,
                                       0.));
  }

  ASSERT_CALL(
    sleqp_vec_create_empty(&linear_lsq_current, linear_lsq_num_variables));

  ASSERT_CALL(
    sleqp_vec_create_empty(&linear_lsq_forward, linear_lsq_num_residuals));

  ASSERT_CALL(sleqp_settings_create(&linear_lsq_settings));

  SleqpLSQCallbacks callbacks = {.set_value       = func_set,
                                 .lsq_residuals   = lsq_residuals,
                                 .lsq_jac_forward = lsq_jac_forward,
                                 .lsq_jac_adjoint = lsq_jac_adjoint,
                                 .cons_val        = NULL,
                                 .cons_jac        = NULL,
                                 .func_free       = NULL};

  ASSERT_CALL(sleqp_lsq_func_create(&linear_lsq_func,
                                    &callbacks,
                                    linear_lsq_num_variables,
                                    linear_lsq_num_constraints,
                                    linear_lsq_num_residuals,
                                    0.,
                                    linear_lsq_settings,
                                    NULL));
}

void
linear_lsq_teardown()
{
  ASSERT_CALL(sleqp_func_release(&linear_lsq_func));

  ASSERT_CALL(sleqp_settings_release(&linear_lsq_settings));

  ASSERT_CALL(sleqp_vec_free(&linear_lsq_current));

  ASSERT_CALL(sleqp_vec_free(&linear_lsq_forward));

  ASSERT_CALL(sleqp_vec_free(&linear_lsq_optimal));
  ASSERT_CALL(sleqp_vec_free(&linear_lsq_initial));

  ASSERT_CALL(sleqp_vec_free(&linear_lsq_cons_ub));
  ASSERT_CALL(sleqp_vec_free(&linear_lsq_cons_lb));

  ASSERT_CALL(sleqp_vec_free(&linear_lsq_var_ub));
  ASSERT_CALL(sleqp_vec_free(&linear_lsq_var_lb));

  ASSERT_CALL(sleqp_vec_free(&linear_lsq_rhs));

  ASSERT_CALL(sleqp_mat_release(&linear_lsq_matrix));
}
