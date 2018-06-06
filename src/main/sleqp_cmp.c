#include "sleqp_cmp.h"

#define ABS(x) (((x) >= 0) ? (x) : (-(x)))

const double eps = 1e-8;

double sleqp_infinity()
{
  return 1e20;
}

SLEQP_Bool sleqp_eq(double x,
                    double y)
{
  return ABS(x - y) <= eps;
}

SLEQP_Bool sleqp_lt(double x,
                    double y)
{
  return ((x) - (y)) < -eps;
}

SLEQP_Bool sleqp_gt(double x,
                    double y)
{
  return ((x) - (y)) > eps;
}

SLEQP_Bool sleqp_le(double x,
                    double y)
{
  return ((x) - (y)) <= eps;
}

SLEQP_Bool sleqp_ge(double x,
                    double y)
{
  return ((x) - (y)) >= -eps;
}

SLEQP_Bool sleqp_neg(double x)
{
  return x < -eps;
}

SLEQP_Bool sleqp_pos(double x)
{
  return x > eps;
}

SLEQP_Bool sleqp_zero(double x)
{
  return ABS(x) <= eps;
}
