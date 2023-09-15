#ifndef SLEQP_MEASURE_H
#define SLEQP_MEASURE_H

#include "direction.h"
#include "iterate.h"
#include "problem.h"
#include "pub_iterate.h"
#include "sparse/pub_vec.h"

#ifdef __cplusplus
extern "C"
{
#endif

  typedef struct SleqpMeasure SleqpMeasure;

  SLEQP_WARNUNUSED
  SLEQP_RETCODE
  sleqp_measure_create(SleqpMeasure** star,
                       SleqpProblem* problem,
                       SleqpSettings* settings);

  SLEQP_WARNUNUSED
  SLEQP_RETCODE
  sleqp_measure_set_iterates(SleqpMeasure* measure,
                             SleqpIterate* iterate,
                             SleqpIterate* trial_iterate,
                             SleqpDirection* direction);

  SLEQP_RETCODE
  sleqp_measure_set_penalty_parameter(SleqpMeasure* measure,
                                      double penalty_parameter);

  SleqpVec*
  sleqp_measure_cons_nonlin(SleqpMeasure* measure);

  double
  sleqp_measure_obj_nonlin(SleqpMeasure* measure);

  double
  sleqp_measure_step_norm(SleqpMeasure* measure);

  SLEQP_WARNUNUSED
  SLEQP_RETCODE
  sleqp_measure_report_trial_point(SleqpMeasure* measure,
                                   const SleqpVec* multipliers);

  SLEQP_WARNUNUSED
  SLEQP_RETCODE
  sleqp_measure_report_soc_trial_point(SleqpMeasure* measure,
                                       SleqpIterate* soc_iterate);

  SLEQP_WARNUNUSED
  SLEQP_RETCODE
  sleqp_measure_lag_nonlin(SleqpMeasure* measure,
                           const SleqpVec* multipliers,
                           double* lag_nonlinearity);

  SLEQP_WARNUNUSED
  SLEQP_RETCODE
  sleqp_measure_capture(SleqpMeasure* measure);

  SLEQP_WARNUNUSED
  SLEQP_RETCODE
  sleqp_measure_release(SleqpMeasure** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_MEASURE_H */
