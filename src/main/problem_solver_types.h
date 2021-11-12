#ifndef PROBLEM_SOLVER_TYPES_H
#define PROBLEM_SOLVER_TYPES_H

#include "iterate.h"

typedef enum {
        SLEQP_PROBLEM_SOLVER_STATUS_UNKNOWN,
        SLEQP_PROBLEM_SOLVER_STATUS_RUNNING,
        SLEQP_PROBLEM_SOLVER_STATUS_OPTIMAL,
        SLEQP_PROBLEM_SOLVER_STATUS_UNBOUNDED,
        SLEQP_PROBLEM_SOLVER_STATUS_LOCALLY_INFEASIBLE,
        SLEQP_PROBLEM_SOLVER_STATUS_ABORT_DEADPOINT,
        SLEQP_PROBLEM_SOLVER_STATUS_ABORT_ITER,
        SLEQP_PROBLEM_SOLVER_STATUS_ABORT_MANUAL,
        SLEQP_PROBLEM_SOLVER_STATUS_ABORT_TIME
} SLEQP_PROBLEM_SOLVER_STATUS;

typedef enum {
        SLEQP_PROBLEM_SOLVER_EVENT_ACCEPTED_ITERATE = 0,
        SLEQP_PROBLEM_SOLVER_EVENT_PERFORMED_ITERATION,
        SLEQP_PROBLEM_SOLVER_NUM_EVENTS
} SLEQP_PROBLEM_SOLVER_EVENT;

typedef struct SleqpProblemSolver SleqpProblemSolver;

typedef SLEQP_RETCODE (*SLEQP_PROBLEM_SOLVER_ACCEPTED_ITERATE)(SleqpProblemSolver* solver,
                                                               SleqpIterate* iterate,
                                                               SleqpIterate* trial_iterate,
                                                               void* callback_data);

typedef SLEQP_RETCODE (*SLEQP_PROBLEM_SOLVER_PERFORMED_ITERATION)(SleqpProblemSolver* solver,
                                                                  void* callback_data);

#endif /* PROBLEM_SOLVER_TYPES_H */
