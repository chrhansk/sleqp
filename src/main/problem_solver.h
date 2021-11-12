#ifndef SLEQP_PROBLEM_SOLVER_H
#define SLEQP_PROBLEM_SOLVER_H

#include "iterate.h"
#include "options.h"

#include "callback_handler.h"
#include "deriv_check.h"
#include "problem_solver_types.h"
#include "trial_point.h"

#include "step/step_rule.h"

struct SleqpProblemSolver
{
        int refcount;

        SleqpProblem* problem;
        SleqpParams* params;
        SleqpOptions* options;

        SLEQP_SOLVER_PHASE solver_phase;

        double* dense_cache;

        SleqpSparseVec* primal_diff;
        SleqpSparseVec* cons_dual_diff;
        SleqpSparseVec* vars_dual_diff;

        SleqpIterate* iterate;
        SleqpIterate* trial_iterate;

        SleqpTimer* elapsed_timer;

        SleqpTrialPointSolver* trial_point_solver;

        SleqpStepRule* step_rule;

        SleqpDerivCheckData* deriv_check;

        SleqpMerit* merit;

        SleqpCallbackHandler* callback_handlers[SLEQP_PROBLEM_SOLVER_NUM_EVENTS];

        SLEQP_PROBLEM_SOLVER_STATUS status;

        SLEQP_STEPTYPE last_step_type;

        double slackness_residuum;

        double stationarity_residuum;

        double feasibility_residuum;

        double trust_radius;

        double lp_trust_radius;

        double penalty_parameter;

        int iteration;

        int elapsed_iterations;

        double time_limit;

        double primal_diff_norm;

        double dual_diff_norm;

        int boundary_step;

        bool abort_next;

        double current_merit_value;

        bool abort_on_local_infeasibility;
};

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_problem_solver_create(SleqpProblemSolver** star,
                                          SLEQP_SOLVER_PHASE solver_phase,
                                          SleqpProblem* problem,
                                          SleqpParams* params,
                                          SleqpOptions* options,
                                          SleqpSparseVec* primal);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_problem_solver_set_primal(SleqpProblemSolver* solver,
                                              const SleqpSparseVec* primal);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_problem_solver_add_callback(SleqpProblemSolver* solver,
                                                SLEQP_PROBLEM_SOLVER_EVENT solver_event,
                                                void* callback_func,
                                                void* callback_data);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_problem_solver_remove_callback(SleqpProblemSolver* solver,
                                                   SLEQP_PROBLEM_SOLVER_EVENT solver_event,
                                                   void* callback_func,
                                                   void* callback_data);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_problem_solver_reset(SleqpProblemSolver* solver);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_problem_solver_abort(SleqpProblemSolver* solver);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_problem_solver_set_iteration(SleqpProblemSolver* solver,
                                                 int iteration);

int sleqp_problem_solver_elapsed_iterations(const SleqpProblemSolver* solver);

double sleqp_problem_solver_get_elapsed_seconds(const SleqpProblemSolver* solver);

SleqpIterate* sleqp_problem_solver_get_iterate(const SleqpProblemSolver* solver);

SLEQP_PROBLEM_SOLVER_STATUS sleqp_problem_solver_get_status(const SleqpProblemSolver* solver);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_problem_solver_print_stats(const SleqpProblemSolver* solver);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_problem_solver_set_func_value(SleqpProblemSolver* solver,
                                                  SleqpIterate* iterate,
                                                  SLEQP_VALUE_REASON reason,
                                                  bool* reject);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_problem_solver_update_trust_radius(SleqpProblemSolver* solver,
                                                       double reduction_ratio,
                                                       bool trial_step_accepted,
                                                       double direction_norm);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_problem_solver_update_lp_trust_radius(SleqpProblemSolver* solver,
                                                          bool trial_step_accepted,
                                                          double trial_step_infnorm,
                                                          double cauchy_step_infnorm,
                                                          bool full_cauchy_step,
                                                          double eps,
                                                          double* lp_trust_radius);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_problem_solver_get_real_state(const SleqpProblemSolver* solver,
                                                  SLEQP_SOLVER_STATE_REAL state,
                                                  double* value);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_problem_solver_get_int_state(const SleqpProblemSolver* solver,
                                                 SLEQP_SOLVER_STATE_INT state,
                                                 int* value);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_problem_solver_get_vec_state(const SleqpProblemSolver* solver,
                                                 SLEQP_SOLVER_STATE_VEC value,
                                                 SleqpSparseVec* result);

SLEQP_RETCODE sleqp_problem_solver_print_header(SleqpProblemSolver* solver);

SLEQP_RETCODE sleqp_problem_solver_print_initial_line(SleqpProblemSolver* solver);

SLEQP_RETCODE sleqp_problem_solver_print_line(SleqpProblemSolver* solver);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_problem_solver_perform_iteration(SleqpProblemSolver* solver);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_problem_solver_solve(SleqpProblemSolver* solver,
                                         int max_num_iterations,
                                         double time_limit,
                                         bool abort_on_local_infeasibility);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_problem_solver_reject_step(SleqpProblemSolver* solver);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_problem_solver_accept_step(SleqpProblemSolver* solver);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_problem_solver_capture(SleqpProblemSolver* solver);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_problem_solver_release(SleqpProblemSolver** star);



#endif /* SLEQP_PROBLEM_SOLVER_H */
