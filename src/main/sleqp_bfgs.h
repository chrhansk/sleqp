#ifndef SLEQP_BFGS_H
#define SLEQP_BFGS_H

/**
 * @file sleqp_bfgs.h
 * @brief Defintion of BFGS method.
 **/

#include "sleqp_options.h"
#include "sleqp_params.h"
#include "sleqp_problem.h"

#include "sleqp_iterate.h"
#include "sparse/sleqp_sparse_vec.h"

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

  typedef struct SleqpBFGSData SleqpBFGSData;

  SLEQP_RETCODE sleqp_bfgs_data_create(SleqpBFGSData** star,
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
  SLEQP_RETCODE sleqp_bfgs_data_push(SleqpBFGSData* data,
                                     SleqpIterate* old_iterate,
                                     SleqpIterate* new_iterate,
                                     SleqpSparseVec* multipliers);

  /**
   * Computes the product of the Hessian approximation with the given direction
   *
   * @param      data         The BFGS data
   * @param[in]  direction    The direction
   * @param[out] product      The product
   *
   **/
  SLEQP_RETCODE sleqp_bfgs_data_hess_prod(const SleqpBFGSData* data,
                                          const SleqpSparseVec* direction,
                                          SleqpSparseVec* product);

  SleqpFunc* sleqp_bfgs_get_func(SleqpBFGSData* data);

  SLEQP_RETCODE sleqp_bfgs_data_capture(SleqpBFGSData* data);

  SLEQP_RETCODE sleqp_bfgs_data_release(SleqpBFGSData** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_BFGS_H */
