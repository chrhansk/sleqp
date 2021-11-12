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

#include "quasi_newton_types.h"

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
SLEQP_NODISCARD
SLEQP_RETCODE sleqp_bfgs_create(SleqpQuasiNewton** star,
                                SleqpFunc* func,
                                SleqpParams* params,
                                SleqpOptions* options);


#endif /* SLEQP_BFGS_H */
