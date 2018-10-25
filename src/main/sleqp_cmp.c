#include "sleqp_cmp.h"

#define ABS(x) (((x) >= 0) ? (x) : (-(x)))

const double eps = 1e-8;

double sleqp_infinity()
{
  return 1e100;
}

bool sleqp_is_inf(double value)
{
  return value / 2. >= sleqp_infinity();
}

bool sleqp_eq(double x,
              double y)
{
  return ABS(x - y) <= eps;
}

bool sleqp_lt(double x,
              double y)
{
  return ((x) - (y)) < -eps;
}

bool sleqp_gt(double x,
              double y)
{
  return ((x) - (y)) > eps;
}

bool sleqp_le(double x,
              double y)
{
  return ((x) - (y)) <= eps;
}

bool sleqp_ge(double x,
              double y)
{
  return ((x) - (y)) >= -eps;
}

bool sleqp_neg(double x)
{
  return x < -eps;
}

bool sleqp_pos(double x)
{
  return x > eps;
}

bool sleqp_zero(double x)
{
  return ABS(x) <= eps;
}

bool sleqp_is_infinity(double x)
{
  return x >= sleqp_infinity();
}
