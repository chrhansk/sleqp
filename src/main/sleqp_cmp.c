#include "sleqp_cmp.h"

const double eps = 1e-16;

double rel_diff(double x, double y)
{
  double d = x - y;

  double m = SLEQP_MAX(SLEQP_ABS(x),
                       SLEQP_ABS(y));

  m = SLEQP_MAX(m, 1.);

  return d / m;
}

double sleqp_eps()
{
  return eps;
}

double sleqp_infinity()
{
  return 1e100;
}

bool sleqp_is_inf(double value)
{
  return value / 2. >= sleqp_infinity();
}

bool sleqp_eq(double x,
              double y,
              double eps)
{
  double d = rel_diff(x, y);
  return SLEQP_ABS(d) <= eps;
}

bool sleqp_lt(double x,
              double y,
              double eps)
{
  return rel_diff(x, y) < -eps;
}

bool sleqp_gt(double x,
              double y,
              double eps)
{
  return rel_diff(x, y) > eps;
}

bool sleqp_le(double x,
              double y,
              double eps)
{
  return rel_diff(x, y) <= eps;
}

bool sleqp_ge(double x,
              double y,
              double eps)
{
  return rel_diff(x, y) >= -eps;
}

// note: these versions are compatible with the ones above, e.g.
// sleqp_neg(x) iff sleqp_lt(x, 0.) and so on

bool sleqp_neg(double x, double eps)
{
  return x < -eps;
}

bool sleqp_pos(double x, double eps)
{
  return x > eps;
}

bool sleqp_zero(double x, double eps)
{
  return SLEQP_ABS(x) <= eps;
}

bool sleqp_is_infinity(double x)
{
  return x >= sleqp_infinity();
}
