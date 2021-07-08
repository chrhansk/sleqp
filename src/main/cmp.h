#ifndef SLEQP_CMP_H
#define SLEQP_CMP_H

#include "pub_cmp.h"

#ifdef __cplusplus
extern "C" {
#endif


  bool sleqp_is_infinite(double value);

  bool sleqp_is_finite(double value);

  bool sleqp_is_eq(double x, double y, double eps);

  bool sleqp_is_lt(double x, double y, double eps);
  bool sleqp_is_gt(double x, double y, double eps);

  bool sleqp_is_leq(double x, double y, double eps);
  bool sleqp_is_geq(double x, double y, double eps);

  bool sleqp_is_neg(double x, double eps);
  bool sleqp_is_pos(double x, double eps);

  bool sleqp_is_zero(double x, double eps);


#define SLEQP_MAX(a,b) (((a)>(b))?(a):(b))

#define SLEQP_MIN(a,b) (((a)<(b))?(a):(b))

#define SLEQP_ABS(a) (((a)>0)?(a):(-(a)))

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_CMP_H */
