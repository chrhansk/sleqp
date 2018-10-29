#ifndef SLEQP_CAUCHY_H
#define SLEQP_CAUCHY_H

#include "sleqp.h"

#include "sleqp_iterate.h"

#include "lp/sleqp_lpi.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpCauchyData SleqpCauchyData;

  SLEQP_RETCODE sleqp_cauchy_data_create(SleqpCauchyData** star,
                                         SleqpProblem* problem,
                                         SleqpParams* params,
                                         SleqpLPi* lp_interface);

  SLEQP_RETCODE sleqp_cauchy_data_free(SleqpCauchyData** star);

  SLEQP_RETCODE sleqp_cauchy_set_iterate(SleqpCauchyData* cauchy_data,
                                         SleqpIterate* iterate,
                                         double trust_radius);

  SLEQP_RETCODE sleqp_cauchy_solve(SleqpCauchyData* cauchy_data,
                                   SleqpSparseVec* gradient,
                                   double penalty);

  SLEQP_RETCODE sleqp_cauchy_get_active_set(SleqpCauchyData* cauchy_data,
                                            SleqpIterate* iterate,
                                            double trust_radius);

  SLEQP_RETCODE sleqp_cauchy_get_direction(SleqpCauchyData* cauchy_data,
                                           SleqpIterate* iterate,
                                           SleqpSparseVec* direction);

  /**
   * Gets the total constraint vioation of the solution at the given iterate.
   * The total violation is defined as the sum over the violations of all
   * violated constraints.
   **/
  SLEQP_RETCODE sleqp_cauchy_get_violation(SleqpCauchyData* cauchy_data,
                                           double* violation);

  /**
   * Computes the Cauchy step by (approximately) minimizing the
   * quadratic penalty along a direction. The direction will
   * be returned scaled and then returned.
   *
   * @param[in]      cauchy_data        Cauchy data
   * @param[in]      iterate            The current iterate
   * @param[in]      penalty_parameter  The penalty parameter \f$ v \f$
   * @param[in]      hessian_direction  The product of the Hessian with the initial direction
   * @param[in,out]  direction          The direction
   * @param[out]     step_length        The computed step length
   *
   **/
  SLEQP_RETCODE sleqp_cauchy_compute_step(SleqpCauchyData* cauchy_data,
                                          SleqpIterate* iterate,
                                          double penalty_parameter,
                                          SleqpSparseVec* hessian_product,
                                          SleqpSparseVec* direction,
                                          double* step_length);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_CAUCHY_H */
