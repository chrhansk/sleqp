#include "sleqp_cmp.h"

#define ABS(x) (((x) >= 0) ? (x) : (-(x)))

const double eps = 1e-8;

double sleqp_infinity()
{
  return 1e20;
}

int sleqp_eq(double x,
             double y)
{
  return ABS(x - y) <= eps;
}
