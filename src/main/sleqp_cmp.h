#ifndef SLEQP_CMP_H
#define SLEQP_CMP_H

#include "sleqp_types.h"

#ifdef __cplusplus
extern "C" {
#endif

  double sleqp_infinity();

  bool sleqp_is_inf(double value);

  bool sleqp_eq(double x, double y);

  bool sleqp_lt(double x, double y);
  bool sleqp_gt(double x, double y);

  bool sleqp_le(double x, double y);
  bool sleqp_ge(double x, double y);

  bool sleqp_neg(double x);
  bool sleqp_pos(double x);

  bool sleqp_zero(double x);

  double sleqp_eps();

  bool sleqp_is_infinity(double x);

#define SLEQP_MAX(a,b) (((a)>(b))?(a):(b))

#define SLEQP_MIN(a,b) (((a)<(b))?(a):(b))

#define SLEQP_ABS(a) (((a)>0)?(a):(-(a)))

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_CMP_H */
