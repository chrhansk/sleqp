#ifndef SLEQP_PUB_TYPES_H
#define SLEQP_PUB_TYPES_H

/**
 * @file pub_types.h
 * @brief Definition of basic types.
 **/

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdbool.h>
#include <stddef.h>

#include "pub_log.h"

  typedef enum {
    SLEQP_OKAY = 0,
    SLEQP_NOMEM,
    SLEQP_FAILED_ASSERTION,
    SLEQP_ILLEGAL_ARGUMENT,
    SLEQP_INVALID_DERIV,
    SLEQP_INTERNAL_ERROR,
    SLEQP_MATH_ERROR
  } SLEQP_RETCODE;

#ifdef SLEQP_HAVE_WARN_UNUSED_RESULT
#define SLEQP_NODISCARD __attribute__((warn_unused_result))
#else
#define SLEQP_NODISCARD
#endif

  typedef enum {
    SLEQP_INACTIVE = 0,
    SLEQP_ACTIVE_LOWER = (1 << 1),
    SLEQP_ACTIVE_UPPER = (1 << 2),
    SLEQP_ACTIVE_BOTH = (SLEQP_ACTIVE_LOWER | SLEQP_ACTIVE_UPPER),
  } SLEQP_ACTIVE_STATE;

  typedef enum {
    SLEQP_OPTIMAL,
    SLEQP_FEASIBLE,
    SLEQP_UNBOUNDED,
    SLEQP_INFEASIBLE,
    SLEQP_INVALID
  } SLEQP_STATUS;

  extern const char* const sleqp_retcode_names[];

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
    SLEQP_STEPTYPE_ACCEPTED_SOC,
    SLEQP_STEPTYPE_REJECTED
  } SLEQP_STEPTYPE;

  typedef enum {
    SLEQP_DUAL_ESTIMATION_TYPE_LP,
    SLEQP_DUAL_ESTIMATION_TYPE_LSQ,
  } SLEQP_DUAL_ESTIMATION_TYPE;

  typedef enum {
    SLEQP_TR_SOLVER_TRLIB = 0,
    SLEQP_TR_SOLVER_CG,
    SLEQP_TR_SOLVER_LSQR,
    SLEQP_TR_SOLVER_AUTO
  } SLEQP_TR_SOLVER;

  typedef enum {
    SLEQP_PARAMETRIC_CAUCHY_DISABLED,
    SLEQP_PARAMETRIC_CAUCHY_COARSE,
    SLEQP_PARAMETRIC_CAUCHY_FINE
  } SLEQP_PARAMETRIC_CAUCHY;

  typedef enum {
    SLEQP_LINESEARCH_EXACT,
    SLEQP_LINESEARCH_APPROX,
  } SLEQP_LINESEARCH;

  typedef enum {
    SLEQP_SOLVER_EVENT_ACCEPTED_ITERATE = 0,
    SLEQP_SOLVER_EVENT_PERFORMED_ITERATION,
    SLEQP_SOLVER_EVENT_FINISHED,
    SLEQP_SOLVER_NUM_EVENTS
  } SLEQP_SOLVER_EVENT;

  typedef enum {
    SLEQP_PREPROCESSING_RESULT_SUCCESS,
    SLEQP_PREPROCESSING_RESULT_FAILURE,
    SLEQP_PREPROCESSING_RESULT_INFEASIBLE
  } SLEQP_PREPROCESSING_RESULT;

  typedef enum {
    SLEQP_STEP_RULE_DIRECT,
    SLEQP_STEP_RULE_WINDOW,
    SLEQP_STEP_RULE_MINSTEP
  } SLEQP_STEP_RULE;

  /**None value to be used in place of integer parameters **/
#define SLEQP_NONE (-1)

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_PUB_TYPES_H */
