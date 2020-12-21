#ifndef SLEQP_TYPES_H
#define SLEQP_TYPES_H

/**
 * @file sleqp_types.h
 * @brief Definition of basic types.
 **/

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdbool.h>
#include <stddef.h>

#include "sleqp_log.h"

  typedef enum {
    SLEQP_OKAY = 0,
    SLEQP_NOMEM,
    SLEQP_FAILED_ASSERTION,
    SLEQP_ILLEGAL_ARGUMENT,
    SLEQP_INVALID_DERIV,
    SLEQP_INTERNAL_ERROR,
    SLEQP_MATH_ERROR
  } SLEQP_RETCODE;

  typedef enum {
    SLEQP_INACTIVE,
    SLEQP_ACTIVE_LOWER,
    SLEQP_ACTIVE_UPPER,
    SLEQP_ACTIVE_BOTH
  } SLEQP_ACTIVE_STATE;

  typedef enum {
    SLEQP_OPTIMAL,
    SLEQP_FEASIBLE,
    SLEQP_INFEASIBLE,
    SLEQP_INVALID
  } SLEQP_STATUS;

  extern const char* sleqp_retcode_names[];

#define SLEQP_CALL(x)                                                      \
  do                                                                       \
  {                                                                        \
    SLEQP_RETCODE status;                                                  \
    if( (status = (x)) != SLEQP_OKAY )                                     \
    {                                                                      \
       sleqp_log_error("Error <%d> (%s) in function %s",                   \
                       status,                                             \
                       sleqp_retcode_names[status],                        \
                       __func__);                                          \
      return status;                                                       \
    }                                                                      \
  }                                                                        \
  while(0)

  typedef enum {
    SLEQP_DERIV_CHECK_SKIP  = (1 << 0),
    SLEQP_DERIV_CHECK_FIRST = (1 << 1),
    SLEQP_DERIV_CHECK_SECOND_EXHAUSTIVE = (1 << 2),
    SLEQP_DERIV_CHECK_SECOND_SIMPLE = (1 << 3),
  } SLEQP_DERIV_CHECK;

  typedef enum {
    SLEQP_HESSIAN_EVAL_EXACT,
    SLEQP_HESSIAN_EVAL_SR1,
    SLEQP_HESSIAN_EVAL_SIMPLE_BFGS,
    SLEQP_HESSIAN_EVAL_DAMPED_BFGS,
  } SLEQP_HESSIAN_EVAL;

  typedef enum {
    SLEQP_BFGS_SIZING_NONE = 0,     // No sizing
    SLEQP_BFGS_SIZING_CENTERED_OL,  // Centered Oren–Luenberger
  } SLEQP_BFGS_SIZING;

  typedef enum {
    SLEQP_STEPTYPE_NONE = 0,
    SLEQP_STEPTYPE_ACCEPTED,
    SLEQP_STEPTYPE_ACCEPTED_FULL,
    SLEQP_STEPTYPE_SOC_ACCEPTED,
    SLEQP_STEPTYPE_REJECTED
  } SLEQP_STEPTYPE;

  typedef enum {
    SLEQP_DUAL_ESTIMATION_TYPE_LP,
    SLEQP_DUAL_ESTIMATION_TYPE_LSQ,
  } SLEQP_DUAL_ESTIMATION_TYPE;

  /**None value to be used in place of integer parameters **/
#define SLEQP_NONE (-1)

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_TYPES_H */
