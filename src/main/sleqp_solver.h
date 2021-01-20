#ifndef SLEQP_SOLVER_H
#define SLEQP_SOLVER_H

/**
 * @file sleqp_solver.h
 * @brief Definition of the solver structure.
 **/

#include "sleqp_export.h"
#include "sleqp_iterate.h"
#include "sleqp_options.h"
#include "sleqp_params.h"
#include "sleqp_problem.h"
#include "sleqp_scale.h"

#ifdef __cplusplus
extern "C" {
#endif

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
  SLEQP_EXPORT SLEQP_RETCODE sleqp_solver_create(SleqpSolver** star,
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
  SLEQP_EXPORT SLEQP_RETCODE sleqp_solver_solve(SleqpSolver* solver,
                                                int max_num_iterations,
                                                double time_limit);

  /**
   * Returns the status of the last call to @ref sleqp_solver_solve
   *
   * @param[in]  solver           The solver
   *
   **/
  SLEQP_EXPORT SLEQP_STATUS sleqp_solver_get_status(SleqpSolver* solver);

  /**
   * Returns the current iterate of the solver
   *
   * @param[in]  solver           The solver
   * @param[out] iterate          A pointer to the current iterate
   *
   **/
  SLEQP_EXPORT SLEQP_RETCODE sleqp_solver_get_solution(SleqpSolver* solver,
                                                       SleqpIterate** iterate);

  SLEQP_EXPORT SLEQP_RETCODE sleqp_solver_get_violated_constraints(SleqpSolver* solver,
                                                                   SleqpIterate* iterate,
                                                                   int* violated_constraints,
                                                                   int* num_violated_constraints);

  /**
   * Returns the number of iterations performed during the last call to @ref sleqp_solver_solve
   *
   * @param[in]  solver           The solver
   *
   **/
  SLEQP_EXPORT int sleqp_solver_get_iterations(SleqpSolver* solver);

  /**
   * Returns the number of seconds elapsed during the last call to @ref sleqp_solver_solve
   *
   * @param[in]  solver           The solver
   *
   **/
  SLEQP_EXPORT double sleqp_solver_get_elapsed_seconds(SleqpSolver* solver);

  SLEQP_EXPORT SLEQP_RETCODE sleqp_solver_capture(SleqpSolver* solver);

  SLEQP_EXPORT SLEQP_RETCODE sleqp_solver_release(SleqpSolver** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_SOLVER_H */
