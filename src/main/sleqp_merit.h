#ifndef SLEQP_MERIT_H
#define SLEQP_MERIT_H

/**
 * @file sleqp_merit.h
 * @brief Definition of merit functions.
 **/

#ifdef __cplusplus
extern "C" {
#endif

#include "sleqp_working_set.h"
#include "sleqp_func.h"
#include "sleqp_iterate.h"
#include "sleqp_params.h"


  typedef struct SleqpMeritData SleqpMeritData;

  SLEQP_NODISCARD SLEQP_RETCODE
  sleqp_merit_data_create(SleqpMeritData** star,
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
  SLEQP_NODISCARD SLEQP_RETCODE
  sleqp_merit_func(SleqpMeritData* merit_data,
                   SleqpIterate* iterate,
                   double penalty_parameter,
                   double* merit_value);

  /**
   * Computes the linearized merit value at the given iterate.
   * The linearized merit value at \f$ \overline{x} \f$ w.r.t.
   * a direction \f$ d \f$ is given by
   *
   * \f[
   * \ell_v(\overline{x}, d) := f(x) + \langle \nabla f(\overline{x}), d \rangle
   * + v \left( \sum_{i=1}^{m} \max((c_i{\overline{x}} + \langle \nabla c_i(\overline{x}) , d \rangle) - u_i, 0) \right)
   * + v  \left(\sum_{i=1}^{m} \max(l_i - (c_i{\overline{x}} + \langle \nabla c_i(\overline{x}) , d \rangle), 0) \right)
   * \f]
   *
   *
   * @param[in]  merit_data        Merit data
   * @param[in]  iterate           The current iterate \f$ overline{x} \f$
   * @param[in]  direction         The direction \f$ d \f$
   * @param[in]  penalty_parameter The penalty parameter \f$ v \f$
   * @param[out] merit_value       The linearized merit value
   *
   **/
  SLEQP_NODISCARD SLEQP_RETCODE
  sleqp_merit_linear(SleqpMeritData* merit_data,
                     SleqpIterate* iterate,
                     const SleqpSparseVec* direction,
                     double penalty_parameter,
                     double* merit_value);

  /**
   * Computes the quadratic merit value at the given iterate.
   * The quadratic merit value at \f$ \overline{x} \f$ w.r.t.
   * a direction \f$ d \f$ is given by
   *
   * \f[
   * q_v(\overline{x}, d, \mu) := \ell_v(\overline{x}, d) + 1/2 \langle d, H(\overline{x}, \mu) d \rangle
   * \f]
   *
   * The computation involves the computation of one Hessian product of the underlying function
   *
   * @param[in]  merit_data        Merit data
   * @param[in]  iterate           The current iterate \f$ overline{x} \f$
   * @param[in]  direction         The direction \f$ d \f$
   * @param[in]  penalty_parameter The penalty parameter \f$ v \f$
   * @param[out] merit_value       The quadratic merit value
   *
   **/
  SLEQP_NODISCARD SLEQP_RETCODE
  sleqp_merit_quadratic(SleqpMeritData* merit_data,
                        SleqpIterate* iterate,
                        const double* func_dual,
                        const SleqpSparseVec* direction,
                        const SleqpSparseVec* cons_duals,
                        double penalty_parameter,
                        double* merit_value);

  SLEQP_NODISCARD SLEQP_RETCODE
  sleqp_merit_data_capture(SleqpMeritData* merit_data);

  SLEQP_NODISCARD SLEQP_RETCODE
  sleqp_merit_data_release(SleqpMeritData** star);


#ifdef __cplusplus
}
#endif

#endif /* SLEQP_MERIT_H */
