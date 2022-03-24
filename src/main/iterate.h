#ifndef SLEQP_ITERATE_H
#define SLEQP_ITERATE_H

#include "pub_iterate.h"

#include "params.h"

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_iterate_slackness_residuum(SleqpProblem* problem,
                                 SleqpIterate* iterate,
                                 double* slackness_residuum);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_iterate_vars_slackness_residuals(SleqpProblem* problem,
                                       SleqpIterate* iterate,
                                       SleqpVec* residuals,
                                       double zero_eps);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_iterate_cons_slackness_residuals(SleqpProblem* problem,
                                       SleqpIterate* iterate,
                                       SleqpVec* residuals,
                                       double zero_eps);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_iterate_feasibility_residuum(SleqpProblem* problem,
                                   SleqpIterate* iterate,
                                   double* feasibility_residuum);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_iterate_stationarity_residuals(SleqpProblem* problem,
                                     SleqpIterate* iterate,
                                     double* cache,
                                     SleqpVec* residuals,
                                     double zero_eps);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_iterate_stationarity_residuum(SleqpProblem* problem,
                                    SleqpIterate* iterate,
                                    double* cache,
                                    double* stationarity_residuum);

bool
sleqp_iterate_is_feasible(SleqpIterate* iterate,
                          double feasibility_residuum,
                          double feasibility_tolerance);

bool
sleqp_iterate_is_optimal(SleqpIterate* iterate,
                         SleqpParams* params,
                         double feasibility_residuum,
                         double slackness_residuum,
                         double stationarity_residuum);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_iterate_get_violated_constraints(SleqpProblem* problem,
                                       SleqpIterate* iterate,
                                       int* violated_constraints,
                                       int* num_violated_constraints,
                                       double feas_eps);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_iterate_copy(const SleqpIterate* source, SleqpIterate* target);

#endif /* SLEQP_ITERATE_H */
