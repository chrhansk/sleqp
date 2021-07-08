#ifndef SLEQP_CAUCHY_H
#define SLEQP_CAUCHY_H

/**
 * @file cauchy.h
 * @brief Definition of Cauchy step-related functions.
 **/

#include "iterate.h"
#include "options.h"
#include "params.h"

#include "lp/lpi.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpCauchy SleqpCauchy;

  typedef enum {
    SLEQP_CAUCHY_OBJECTIVE_TYPE_DEFAULT,
    SLEQP_CAUCHY_OBJECTIVE_TYPE_FEASIBILITY,
    SLEQP_CAUCHY_OBJECTIVE_TYPE_MIXED,
    SLEQP_NUM_CAUCHY_OBJECTIVES
  } SLEQP_CAUCHY_OBJECTIVE_TYPE;

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_cauchy_create(SleqpCauchy** star,
                                    SleqpProblem* problem,
                                    SleqpParams* params,
                                    SleqpOptions* options,
                                    SleqpLPi* lp_interface);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_cauchy_release(SleqpCauchy** star);

  /**
   * Sets the iterate and trust radius for the current LP.
   *
   * @param[in]       cauchy_data        Cauchy data
   * @param[in]       iterate            The current iterate
   * @param[in]       trust_radius       The trust radius
   **/
  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_cauchy_set_iterate(SleqpCauchy* cauchy_data,
                                         SleqpIterate* iterate,
                                         double trust_radius);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_cauchy_set_trust_radius(SleqpCauchy* cauchy_data,
                                              double trust_radius);

  /**
   * Solve the LP according to the current iterate. The
   * objective is determined according to the function
   * gradient and penalty value.
   *
   * @param[in]       cauchy_data        Cauchy data
   * @param[in]       gradient           The objective value gradient
   * @param[in]       penalty            The penalty value
   *
   **/
  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_cauchy_solve(SleqpCauchy* cauchy_data,
                                   SleqpSparseVec* gradient,
                                   double penalty,
                                   SLEQP_CAUCHY_OBJECTIVE_TYPE objective_type);

  /**
   * Returns the value of the linearized merit function at the
   * optimal solution of the current LP.
   *
   * @param[in]       cauchy_data        Cauchy data
   * @param[in]       objective_value    The objective value
   *
   **/
  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_cauchy_get_objective_value(SleqpCauchy* cauchy_data,
                                                 double* objective_value);

  /**
   * Updates the working of the given iterate according to the optimal basis of
   * the current LP.
   *
   * @param[in]       cauchy_data        Cauchy data
   * @param[out]      iterate            The current iterate
   *
   **/
  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_cauchy_get_working_set(SleqpCauchy* cauchy_data,
                                             SleqpIterate* iterate);

  /**
   * Returns the Cauchy direction according to the current LP solution.
   *
   * @param[in]       cauchy_data        Cauchy data
   * @param[out]      direction          The Cauchy direction
   *
   **/
  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_cauchy_get_direction(SleqpCauchy* cauchy_data,
                                           SleqpSparseVec* direction);

  /**
   * Determines whether the current iterate is locally infeasible, i.e.,
   * no direction to feasibility in the linearization can be found while
   * none of the trust region constraints are active.
   *
   * @param[in]       cauchy_data          Cauchy data
   * @param[in]       iterate              The current iterate
   * @param[in]       trust_radius         The trust radius
   * @param[out]      locally_infeasible   Whether or not the iterate is locally infeasible
   **/
  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_cauchy_locally_infeasible(SleqpCauchy* cauchy_data,
                                                bool* locally_infeasible);

  /**
   * Returns the dual estimation according to duals of the current LP solution. Duals
   * will be stored in the given iterate.
   *
   * @param[in]       cauchy_data        Cauchy data
   * @param[out]      iterate            The current iterate
   *
   **/
  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_cauchy_get_dual_estimation(SleqpCauchy* cauchy_data,
                                                 SleqpIterate* iterate);

  /**
   * Gets the total constraint vioation of the solution at the given iterate.
   * The total violation is defined as the sum over the violations of all
   * violated constraints.
   *
   * @param[in]       cauchy_data        Cauchy data
   * @param[out]      violation          The total constraint violation
   **/
  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_cauchy_get_violation(SleqpCauchy* cauchy_data,
                                           double* violation);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_CAUCHY_H */
