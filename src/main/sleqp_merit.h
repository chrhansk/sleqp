#ifndef SLEQP_MERIT_H
#define SLEQP_MERIT_H

/**
 * @file sleqp_merit.h
 * @brief Definition of merit functions.
 **/

#ifdef __cplusplus
extern "C" {
#endif

#include "sleqp_active_set.h"
#include "sleqp_func.h"
#include "sleqp_iterate.h"
#include "sleqp_params.h"
#include "sleqp.h"


  typedef struct SleqpMeritData SleqpMeritData;

  SLEQP_RETCODE sleqp_merit_data_create(SleqpMeritData** star,
                                        SleqpProblem* problem,
                                        SleqpParams* params);

  /**
   * Computes the exact merit value at the given
   * iterate. The exact merit value at \f$ x \f$ is given by
   *
   * \f[
   * \Phi_v(x) := f(x)
   *              + v \left( \sum_{i=1}^{m} \max(c_i(x) - u_i, 0) \right)
   *              + v \left(\sum_{i=1}^{m} \max(l_i - c_i(x), 0) \right)
   * \f]
   *
   * @param[in]  merit_data        Merit data
   * @param[in]  iterate           The current iterate
   * @param[in]  penalty_parameter The penalty parameter \f$ v \f$
   * @param[out] merit_value       The exact merit value
   *
   **/
  SLEQP_RETCODE sleqp_merit_func(SleqpMeritData* merit_data,
                                 SleqpIterate* iterate,
                                 double penalty_parameter,
                                 double* merit_value);

  /**
   * Computes the linearized merit value at the given iterate.
   * The linearized merit value at \f$ \overline{x} \f$ w.r.t.
   * a direction \f$ d \f$ is given by
   *
   * \f[
   * \ell_v(\overline{x}, d) := \langle \nabla f(\overline{x}), d \rangle
   * + v \left( \sum_{i=1}^{m} \max((c_i{\overline{x}} + \langle \nabla c_i(\overline{x}) , d \rangle) - u_i, 0) \right)
   * + v  \left(\sum_{i=1}^{m} \max(l_i - (c_i{\overline{x}} + \langle \nabla c_i(\overline{x}) , d \rangle), 0) \right)
   * \f]
   *
   *
   * @param[in]  merit_data        Merit data
   * @param[in]  iterate           The current iterate \f$ overline{x} \f$
   * @param[in]  direction         The direction \f$ d \f$
   * @param[in]  penalty_parameter The penalty parameter \f$ v \f$
   * @param[out] merit_value       The linear merit value
   *
   **/
  SLEQP_RETCODE sleqp_merit_linear(SleqpMeritData* merit_data,
                                   SleqpIterate* iterate,
                                   SleqpSparseVec* direction,
                                   double penalty_parameter,
                                   double* merit_value);

  SLEQP_RETCODE sleqp_merit_linear_gradient(SleqpMeritData* merit_data,
                                            SleqpIterate* iterate,
                                            SleqpSparseVec* direction,
                                            double penalty_parameter,
                                            SleqpSparseVec* gradient);

  SLEQP_RETCODE sleqp_merit_data_free(SleqpMeritData** star);


#ifdef __cplusplus
}
#endif

#endif /* SLEQP_MERIT_H */
