#include "cmp.h"

#include <math.h>

const double inf = 1e100;

static double rel_diff(double x, double y)
{
  double d = x - y;

  double m = SLEQP_MAX(SLEQP_ABS(x),
                       SLEQP_ABS(y));

  m = SLEQP_MAX(m, 1.);

  return d / m;
}

double sleqp_infinity()
{
  return inf;
}

bool sleqp_is_infinite(double value)
{
  if(isnan(value))
  {
    return false;
  }

  return value >= sleqp_infinity() / 2.;
}

bool sleqp_is_finite(double value)
{
  return !(sleqp_is_infinite(SLEQP_ABS(value)));
}

bool sleqp_is_eq(double x,
                 double y,
                 double eps)
{
  double d = rel_diff(x, y);
  return SLEQP_ABS(d) <= eps;
}

bool sleqp_is_lt(double x,
                 double y,
                 double eps)
{
  return rel_diff(x, y) < -eps;
}

bool sleqp_is_gt(double x,
                 double y,
                 double eps)
{
  return rel_diff(x, y) > eps;
}

bool sleqp_is_leq(double x,
                  double y,
                  double eps)
{
  return rel_diff(x, y) <= eps;
}

bool sleqp_is_geq(double x,
                  double y,
                  double eps)
{
  return rel_diff(x, y) >= -eps;
}

// note: these versions are compatible with the ones above, e.g.
// sleqp_is_neg(x) iff sleqp_is_lt(x, 0.) and so on

bool sleqp_is_neg(double x, double eps)
{
  return x < -eps;
}

bool sleqp_is_pos(double x, double eps)
{
  return x > eps;
}

bool sleqp_is_zero(double x, double eps)
{
  return SLEQP_ABS(x) <= eps;
}
