#ifndef SLEQP_CMP_H
#define SLEQP_CMP_H

/**
 * @file sleqp_cmp.h
 * @brief Definition of numerical comparison functions.
 **/

#include "sleqp_types.h"

#ifdef __cplusplus
extern "C" {
#endif

  double sleqp_infinity();

  bool sleqp_is_inf(double value);

  bool sleqp_eq(double x, double y, double eps);

  bool sleqp_lt(double x, double y, double eps);
  bool sleqp_gt(double x, double y, double eps);

  bool sleqp_le(double x, double y, double eps);
  bool sleqp_ge(double x, double y, double eps);

  bool sleqp_neg(double x, double eps);
  bool sleqp_pos(double x, double eps);

  bool sleqp_zero(double x, double eps);

  double sleqp_eps();


#define SLEQP_MAX(a,b) (((a)>(b))?(a):(b))

#define SLEQP_MIN(a,b) (((a)<(b))?(a):(b))

#define SLEQP_ABS(a) (((a)>0)?(a):(-(a)))

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_CMP_H */
