#ifndef SLEQP_SOLVER_H
#define SLEQP_SOLVER_H

/**
 * @file sleqp_solver.h
 * @brief Definition of the solver structure.
 **/

#include "sleqp_iterate.h"
#include "sleqp_options.h"
#include "sleqp_params.h"
#include "sleqp_problem.h"
#include "sleqp_scale.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpSolver SleqpSolver;

  SLEQP_RETCODE sleqp_solver_create(SleqpSolver** star,
                                    SleqpProblem* problem,
                                    SleqpParams* params,
                                    SleqpOptions* options,
                                    SleqpSparseVec* x,
                                    SleqpScalingData* scaling_data);

  SLEQP_RETCODE sleqp_solver_solve(SleqpSolver* solver,
                                   int max_num_iterations,
                                   double time_limit);

  SLEQP_STATUS sleqp_solver_get_status(SleqpSolver* solver);

  SLEQP_RETCODE sleqp_solver_get_solution(SleqpSolver* solver,
                                          SleqpIterate** iterate);

  SLEQP_RETCODE sleqp_solver_get_violated_constraints(SleqpSolver* solver,
                                                      SleqpIterate* iterate,
                                                      int* violated_constraints,
                                                      int* num_violated_constraints);

  int sleqp_solver_get_iterations(SleqpSolver* solver);

  double sleqp_solver_get_elapsed_seconds(SleqpSolver* solver);

  SLEQP_RETCODE sleqp_solver_free(SleqpSolver** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_SOLVER_H */
