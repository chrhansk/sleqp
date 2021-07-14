#ifndef SLEQP_BFGS_H
#define SLEQP_BFGS_H

/**
 * @file bfgs.h
 * @brief Defintion of BFGS method.
 **/

#include "iterate.h"
#include "options.h"
#include "params.h"
#include "problem.h"
#include "timer.h"

#include "sparse/sparse_vec.h"

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * @defgroup BFGS method
   *
   * The BFGS method is a quasi-Newton method based on
   * a sequence of vectors \f$ (s_k, y_k) \f$, where
   * the former are differences between primal points
   * of iterated and the latter differences in the
   * corresponding Lagrangean gradients
   * \f$ \nabla_x L(x, \lambda, \mu) \f$.
   *
   **/

  typedef struct SleqpBFGS SleqpBFGS;

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_bfgs_create(SleqpBFGS** star,
                                  SleqpFunc* func,
                                  SleqpParams* params,
                                  SleqpOptions* options);

  /**
   * Pushes a new pair \f$ (s_k, y_k) \f$. If the new pair
   * increases the number of pairs beyond the limit
   * (@see sleqp_options_get_quasi_newton_num_iterates)
   * the oldest pair will be removed. The vecto
   * \f$ y_k \f$ is determined according to the
   * prodided multipliers.
   *
   * @param data         The BFGS data
   * @param old_iterate  The previous iterate
   * @param new_iterate  The current iterate
   * @param multipliers  The multipliers used to compute \f$ y_k \f$
   *
   **/
  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_bfgs_push(SleqpBFGS* data,
                                SleqpIterate* old_iterate,
                                SleqpIterate* new_iterate,
                                SleqpSparseVec* multipliers);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_bfgs_reset(SleqpBFGS* data);

  /**
   * Computes the product of the Hessian approximation with the given direction
   *
   * @param      data         The BFGS data
   * @param[in]  direction    The direction
   * @param[out] product      The product
   *
   **/
  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_bfgs_hess_prod(const SleqpBFGS* data,
                                     const SleqpSparseVec* direction,
                                     SleqpSparseVec* product);

  SleqpFunc* sleqp_bfgs_get_func(SleqpBFGS* data);

  SleqpTimer* sleqp_bfgs_update_timer(SleqpBFGS* data);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_bfgs_capture(SleqpBFGS* data);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_bfgs_release(SleqpBFGS** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_BFGS_H */
