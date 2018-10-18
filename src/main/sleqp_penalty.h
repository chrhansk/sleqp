#ifndef SLEQP_PENALTY_H
#define SLEQP_PENALTY_H

/**
 * @file sleqp_penalty.h
 * @brief Definition of penalty functions.
 **/

#ifdef __cplusplus
extern "C" {
#endif

#include "sleqp_active_set.h"
#include "sleqp_func.h"
#include "sleqp_iterate.h"
#include "sleqp.h"


  typedef struct SleqpPenalty SleqpPenalty;

  SLEQP_RETCODE sleqp_penalty_create(SleqpPenalty** star,
                                     SleqpProblem* problem,
                                     SleqpFunc* func);

  /**
   * Computes the exact penalty value at the given
   * iterate. The exact penalty value at \f$ x \f$ is given by
   *
   * \f[
   * \Phi_v(x) := f(x)
   *              + v \left( \sum_{i=1}^{m} \max(c_i(x) - u_i, 0) \right)
   *              + v \left(\sum_{i=1}^{m} \max(l_i - c_i(x), 0) \right)
   * \f]
   *
   * @param[in]  penalty_data      Penalty data
   * @param[in]  iterate           The current iterate
   * @param[in]  penalty_parameter The penalty parameter \f$ v \f$
   * @param[out] penalty_value     The exact penalty value
   *
   **/
  SLEQP_RETCODE sleqp_penalty_func(SleqpPenalty* penalty_data,
                                   SleqpIterate* iterate,
                                   double penalty_parameter,
                                   double* penalty_value);

  /**
   * Computes the linearized penalty value at the given iterate.
   * The linearized penalty value at \f$ \overline{x} \f$ w.r.t.
   * a direction \f$ d \f$ is given by
   *
   * \f[
   * \ell_v(\overline{x}, d) := \langle \nabla f(\overline{x}), d \rangle
   * + v \left( \sum_{i=1}^{m} \max((c_i{\overline{x}} + \langle \nabla c_i(\overline{x}) , d \rangle) - u_i, 0) \right)
   * + v  \left(\sum_{i=1}^{m} \max(l_i - (c_i{\overline{x}} + \langle \nabla c_i(\overline{x}) , d \rangle), 0) \right)
   * \f]
   *
   *
   * @param[in]  penalty_data      Penalty data
   * @param[in]  iterate           The current iterate \f$ overline{x} \f$
   * @param[in]  direction         The direction \f$ d \f$
   * @param[in]  penalty_parameter The penalty parameter \f$ v \f$
   * @param[out] penalty_value     The linear penalty value
   *
   **/
  SLEQP_RETCODE sleqp_penalty_linear(SleqpPenalty* penalty_data,
                                     SleqpIterate* iterate,
                                     SleqpSparseVec* direction,
                                     double penalty_parameter,
                                     double* penalty_value);

  SLEQP_RETCODE sleqp_penalty_free(SleqpPenalty** star);


#ifdef __cplusplus
}
#endif

#endif /* SLEQP_PENALTY_H */
