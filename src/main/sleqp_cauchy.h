#ifndef SLEQP_CAUCHY_H
#define SLEQP_CAUCHY_H

#include "sleqp.h"

#include "sleqp_active_set.h"
#include "sleqp_iterate.h"

#include "lp/sleqp_lpi.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpCauchyData SleqpCauchyData;

  SLEQP_RETCODE sleqp_cauchy_data_create(SleqpCauchyData** star,
                                         SleqpProblem* problem,
                                         SleqpLPi* lp_interface);

  SLEQP_RETCODE sleqp_cauchy_data_free(SleqpCauchyData** star);

  SLEQP_RETCODE sleqp_cauchy_compute_direction(SleqpCauchyData* cauchy_data,
                                               SleqpIterate* iterate,
                                               double penalty,
                                               double trust_radius);

  SLEQP_RETCODE sleqp_cauchy_get_active_set(SleqpCauchyData* cauchy_data,
                                            SleqpIterate* iterate,
                                            SleqpActiveSet* active_set,
                                            double trust_radius);

  SLEQP_RETCODE sleqp_cauchy_get_direction(SleqpCauchyData* cauchy_data,
                                           SleqpIterate* iterate,
                                           SleqpSparseVec* direction);

  /**
   * Computes the Cauchy step by (approximately) minimizing the
   * quadratic penalty along a direction. The direction will
   * be returned scaled and then returned.
   *
   * @param[in]      cauchy_data        Cauchy data
   * @param[in]      iterate            The current iterate
   * @param[in]      active_set         The active set at the current iterate
   * @param[in]      penalty_parameter  The penalty parameter \f$ v \f$
   * @param[in,out]  direction          The direction
   *
   **/
  SLEQP_RETCODE sleqp_cauchy_compute_step(SleqpCauchyData* cauchy_data,
                                          SleqpIterate* iterate,
                                          SleqpActiveSet* active_set,
                                          double penalty_parameter,
                                          SleqpSparseVec* direction,
                                          double* predicted_reduction);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_CAUCHY_H */
