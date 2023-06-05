#include "tr_util.h"

#include <math.h>

#include "cmp.h"
#include "fail.h"

SLEQP_RETCODE
sleqp_tr_compute_bdry_sol(const SleqpVec* previous,
                          const SleqpVec* direction,
                          SleqpSettings* settings,
                          double radius,
                          SleqpVec* result)
{
  const double zero_eps = sleqp_settings_real_value(settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  const double eps = sleqp_settings_real_value(settings, SLEQP_SETTINGS_REAL_EPS);

  SLEQP_NUM_ASSERT_PARAM(eps);

  double prev_dot_d = 0.;

  SLEQP_CALL(sleqp_vec_dot(previous, direction, &prev_dot_d));

  const double p_norm = sleqp_vec_norm(previous);
  const double d_norm = sleqp_vec_norm(direction);

  const double d_norm_sq = d_norm * d_norm;
  const double p_norm_sq = p_norm * p_norm;

  sleqp_assert_is_leq(p_norm, radius, eps);

  assert(d_norm_sq > 0.);

  const double squared_rad = radius * radius;

  const double inner
    = (prev_dot_d * prev_dot_d) - d_norm_sq * (p_norm_sq - squared_rad);

  assert(inner >= 0.);

  const double factor = 1. / d_norm_sq * (-prev_dot_d + sqrt(inner));

  SLEQP_CALL(
    sleqp_vec_add_scaled(previous, direction, 1., factor, zero_eps, result));

  sleqp_num_assert(sleqp_is_eq(sleqp_vec_norm(result), radius, eps));

  return SLEQP_OKAY;
}
