#ifndef SLEQP_SOLVER_H
#define SLEQP_SOLVER_H

/**
 * @file sleqp_solver.h
 * @brief Definition of the solver structure.
 **/

#include "sleqp_callback_types.h"
#include "sleqp_export.h"
#include "sleqp_iterate.h"
#include "sleqp_options.h"
#include "sleqp_params.h"
#include "sleqp_problem.h"
#include "sleqp_scale.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef enum {
    SLEQP_SOLVER_STATE_REAL_TRUST_RADIUS,
    SLEQP_SOLVER_STATE_REAL_LP_TRUST_RADIUS,
    SLEQP_SOLVER_STATE_REAL_SCALED_FUNC_VAL,
    SLEQP_SOLVER_STATE_REAL_SCALED_MERIT_VAL,
    SLEQP_SOLVER_STATE_REAL_SCALED_FEAS_RES,
    SLEQP_SOLVER_STATE_REAL_SCALED_STAT_RES,
    SLEQP_SOLVER_STATE_REAL_SCALED_SLACK_RES,
    SLEQP_SOLVER_STATE_REAL_PENALTY_PARAM,
  } SLEQP_SOLVER_STATE_REAL;

  typedef enum {
    SLEQP_SOLVER_STATE_INT_LAST_STEP_ON_BDRY,
    SLEQP_SOLVER_STATE_INT_LAST_STEP_TYPE,
  } SLEQP_SOLVER_STATE_INT;

  typedef enum {
    SLEQP_SOLVER_STATE_VEC_SCALED_STAT_RESIDUALS,
    SLEQP_SOLVER_STATE_VEC_SCALED_FEAS_RESIDUALS,
    SLEQP_SOLVER_STATE_VEC_SCALED_CONS_SLACK_RESIDUALS,
    SLEQP_SOLVER_STATE_VEC_SCALED_VAR_SLACK_RESIDUALS,
  } SLEQP_SOLVER_STATE_VEC;

  typedef struct SleqpSolver SleqpSolver;

  /**
   * Creates a solver
   *
   * @param[out] star            A pointer to the solver to be created
   * @param[in]  problem         The underlying problem
   * @param[in]  params          The problem parameters
   * @param[in]  options         The solver options
   * @param[in]  x               The initial solution
   * @param[in]  scaling_data    The scaling to be used (may be `NULL`)
   **/
  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_solver_create(SleqpSolver** star,
                                    SleqpProblem* problem,
                                    SleqpParams* params,
                                    SleqpOptions* options,
                                    SleqpSparseVec* x,
                                    SleqpScalingData* scaling_data);

  /**
   * Solves the problem by performing iteration starting from the current solution
   *
   * @param[in]  solver           The solver
   * @param[in]  num_iterations   The number of iterations to be performed, or @ref SLEQP_NONE
   * @param[in]  time_limit       A time limit in seconds, or @ref SLEQP_NONE
   **/
  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_solver_solve(SleqpSolver* solver,
                                   int max_num_iterations,
                                   double time_limit);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_solver_get_real_state(const SleqpSolver* solver,
                                            SLEQP_SOLVER_STATE_REAL state,
                                            double* value);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_solver_get_int_state(const SleqpSolver* solver,
                                           SLEQP_SOLVER_STATE_INT state,
                                           int* value);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_solver_get_vec_state(const SleqpSolver* solver,
                                           SLEQP_SOLVER_STATE_VEC value,
                                           SleqpSparseVec* result);

  /**
   * Returns the status of the last call to @ref sleqp_solver_solve
   *
   * @param[in]  solver           The solver
   *
   **/
  SLEQP_EXPORT SLEQP_STATUS sleqp_solver_get_status(const SleqpSolver* solver);

  SLEQP_EXPORT SLEQP_RETCODE sleqp_solver_abort(SleqpSolver* solver);

  /**
   * Returns the current iterate of the solver
   *
   * @param[in]  solver           The solver
   * @param[out] iterate          A pointer to the current iterate
   *
   **/
  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_solver_get_solution(SleqpSolver* solver,
                                          SleqpIterate** iterate);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_solver_get_violated_constraints(SleqpSolver* solver,
                                                      SleqpIterate* iterate,
                                                      int* violated_constraints,
                                                      int* num_violated_constraints);

  /**
   * Returns the number of iterations performed during the last call to @ref sleqp_solver_solve
   *
   * @param[in]  solver           The solver
   *
   **/
  SLEQP_EXPORT int sleqp_solver_get_iterations(const SleqpSolver* solver);

  /**
   * Returns the number of seconds elapsed during the last call to @ref sleqp_solver_solve
   *
   * @param[in]  solver           The solver
   *
   **/
  SLEQP_EXPORT double sleqp_solver_get_elapsed_seconds(const SleqpSolver* solver);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_solver_add_callback(SleqpSolver* solver,
                                          SLEQP_SOLVER_EVENT solver_event,
                                          void* callback_func,
                                          void* callback_data);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_solver_remove_callback(SleqpSolver* solver,
                                             SLEQP_SOLVER_EVENT solver_event,
                                             void* callback_func,
                                             void* callback_data);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_solver_capture(SleqpSolver* solver);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_solver_release(SleqpSolver** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_SOLVER_H */
