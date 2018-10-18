#ifndef SLEQP_CMP_H
#define SLEQP_CMP_H

#include "sleqp_types.h"

#ifdef __cplusplus
extern "C" {
#endif

  double sleqp_infinity();

  SLEQP_Bool sleqp_eq(double x, double y);

  SLEQP_Bool sleqp_lt(double x, double y);
  SLEQP_Bool sleqp_gt(double x, double y);

  SLEQP_Bool sleqp_le(double x, double y);
  SLEQP_Bool sleqp_ge(double x, double y);

  SLEQP_Bool sleqp_neg(double x);
  SLEQP_Bool sleqp_pos(double x);

  SLEQP_Bool sleqp_zero(double x);

  SLEQP_Bool sleqp_is_infinity(double x);

#define SLEQP_MAX(a,b) (((a)>(b))?(a):(b))
#define SLEQP_MIN(a,b) (((a)<(b))?(a):(b))

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_CMP_H */
