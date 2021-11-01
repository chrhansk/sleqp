#ifndef SLEQP_SOLVER_H
#define SLEQP_SOLVER_H

#include "pub_solver.h"

#include "callback_handler.h"
#include "deriv_check.h"
#include "merit.h"
#include "polish.h"
#include "problem_scaling.h"
#include "trial_point.h"

#include "quasi_newton/quasi_newton.h"
#include "step/step_rule.h"

#include "preprocessor/preprocessor.h"

#ifdef __cplusplus
extern "C" {
#endif

  struct SleqpSolver
  {
    int refcount;

    SleqpProblem* original_problem;

    SleqpScaling* scaling_data;

    SleqpSparseVec* scaled_primal;
    SleqpSparseVec* primal;

    SleqpProblem* scaled_problem;

    SleqpPreprocessor* preprocessor;

    SleqpProblemScaling* problem_scaling;

    bool restore_original_iterate;
    SleqpIterate* original_iterate;

    SleqpIterate* scaled_iterate;

    SleqpProblem* problem;

    SleqpTimer* elapsed_timer;

    SLEQP_STATUS status;

    SleqpParams* params;

    SleqpOptions* options;

    SleqpDerivCheckData* deriv_check;

    SleqpTrialPointSolver* trial_point_solver;

    SleqpStepRule* step_rule;

    SleqpIterate* iterate;

    SleqpIterate* trial_iterate;

    SleqpSparseVec* original_violation;

    SLEQP_STEPTYPE last_step_type;

    SleqpMerit* merit;

    SleqpPolishing* polishing;

    SleqpCallbackHandler** callback_handlers;

    // Primal / dual step lengths

    SleqpSparseVec* primal_diff;

    double primal_diff_norm;

    SleqpSparseVec* cons_dual_diff;

    SleqpSparseVec* vars_dual_diff;

    double dual_diff_norm;

    double current_merit_value;

    double* dense_cache;

    // residuum

    double slackness_residuum;

    double stationarity_residuum;

    double feasibility_residuum;

    SleqpQuasiNewton* quasi_newton;

    // parameters, adjusted throughout...

    bool locally_infeasible;

    double trust_radius;

    double lp_trust_radius;

    double penalty_parameter;

    // misc

    int boundary_step;

    double elapsed_seconds;

    int iteration;

    double time_limit;

    bool abort_next;
  };

  SLEQP_RETCODE sleqp_solver_print_header(SleqpSolver* solver);

  SLEQP_RETCODE sleqp_solver_print_initial_line(SleqpSolver* solver);

  SLEQP_RETCODE sleqp_solver_print_line(SleqpSolver* solver);

  SLEQP_RETCODE sleqp_solver_print_stats(SleqpSolver* solver,
                                         double violation);

  SLEQP_RETCODE sleqp_solver_update_trust_radius(SleqpSolver* solver,
                                                 double reduction_ratio,
                                                 bool trial_step_accepted,
                                                 double direction_norm);

  SLEQP_RETCODE sleqp_solver_update_lp_trust_radius(SleqpSolver* solver,
                                                    bool trial_step_accepted,
                                                    double trial_step_infnorm,
                                                    double cauchy_step_infnorm,
                                                    bool full_cauchy_step,
                                                    double eps,
                                                    double* lp_trust_radius);

  double sleqp_solver_remaining_time(SleqpSolver* solver);

  SLEQP_RETCODE sleqp_solver_restore_original_iterate(SleqpSolver* solver);

  SLEQP_RETCODE sleqp_solver_perform_iteration(SleqpSolver* solver);

  SLEQP_RETCODE sleqp_solver_set_func_value(SleqpSolver* solver,
                                            SleqpIterate* iterate,
                                            SLEQP_VALUE_REASON reason,
                                            bool* reject);

  SLEQP_RETCODE sleqp_solver_accept_step(SleqpSolver* solver);

  SLEQP_RETCODE sleqp_solver_reject_step(SleqpSolver* solver);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_SOLVER_H */
