#ifndef SLEQP_MERIT_H
#define SLEQP_MERIT_H

/**
 * @file merit.h
 * @brief Definition of merit functions.
 **/

#include "direction.h"
#include "func.h"
#include "iterate.h"
#include "working_set.h"

typedef struct SleqpMeritData SleqpMerit;

SLEQP_NODISCARD SLEQP_RETCODE
sleqp_merit_create(SleqpMerit** star,
                   SleqpProblem* problem,
                   SleqpSettings* params);

/**
 * Computes the exact merit value at the given
 * iterate. The exact merit value at \f$ x \f$ is given by
 *
 * \f[
 * \Phi_{\nu}(x) := f(x)
 *              + {\nu} \left( \sum_{i=1}^{m} \max(c_i(x) - u_i, 0) \right)
 *              + {\nu} \left(\sum_{i=1}^{m} \max(l_i - c_i(x), 0) \right)
 * \f]
 *
 * @param[in]  merit             Merit data
 * @param[in]  iterate           The current iterate
 * @param[in]  penalty_parameter The penalty parameter \f$ \nu \f$
 * @param[out] merit_value       The exact merit value
 *
 **/
SLEQP_NODISCARD SLEQP_RETCODE
sleqp_merit_func(SleqpMerit* merit,
                 SleqpIterate* iterate,
                 double penalty_parameter,
                 double* merit_value);

/**
 * Computes the linearized merit value at the given iterate.
 * The linearized merit value at \f$ \overline{x} \f$ w.r.t.
 * a direction \f$ d \f$ is given by
 *
 * \f[
 * \ell_v(\overline{x}, d) := f(\overline{x}) + \langle \nabla f(\overline{x}), d \rangle
 * + v \left( \sum_{i=1}^{m} \max((c_i(\overline{x}) + \langle \nabla
 * c_i(\overline{x}) , d \rangle) - u_i, 0) \right)
 * + v  \left(\sum_{i=1}^{m} \max(l_i - (c_i(\overline{x}) + \langle \nabla
 * c_i(\overline{x}) , d \rangle), 0) \right)
 * \f]
 *
 *
 * @param[in]  merit             Merit data
 * @param[in]  iterate           The current iterate \f$ overline{x} \f$
 * @param[in]  direction         The direction \f$ d \f$
 * @param[in]  penalty_parameter The penalty parameter \f$ \nu \f$
 * @param[out] merit_value       The linearized merit value
 *
 **/
SLEQP_NODISCARD SLEQP_RETCODE
sleqp_merit_linear(SleqpMerit* merit,
                   SleqpIterate* iterate,
                   SleqpDirection* direction,
                   double penalty_parameter,
                   double* merit_value);

/**
 * Computes the quadratic merit value at the given iterate.
 * The quadratic merit value at \f$ \overline{x} \f$ w.r.t.
 * a direction \f$ d \f$ is given by
 *
 * \f[
 * q_v(\overline{x}, d, \mu) := \ell_v(\overline{x}, d) + 1/2 \langle d,
 * H(\overline{x}, \mu) d \rangle
 * \f]
 *
 * The computation involves the computation of one Hessian product of the
 *underlying function
 *
 * @param[in]  merit             Merit data
 * @param[in]  iterate           The current iterate \f$ overline{x} \f$
 * @param[in]  direction         The direction \f$ d \f$
 * @param[in]  penalty_parameter The penalty parameter \f$ \nu \f$
 * @param[out] merit_value       The quadratic merit value
 *
 **/
SLEQP_NODISCARD SLEQP_RETCODE
sleqp_merit_quadratic(SleqpMerit* merit,
                      SleqpIterate* iterate,
                      SleqpDirection* direction,
                      double penalty_parameter,
                      double* merit_value);

SLEQP_NODISCARD SLEQP_RETCODE
sleqp_merit_capture(SleqpMerit* merit);

SLEQP_NODISCARD SLEQP_RETCODE
sleqp_merit_release(SleqpMerit** star);

#endif /* SLEQP_MERIT_H */
