#ifndef SLEQP_SOLVER_H
#define SLEQP_SOLVER_H

#include "pub_solver.h"

#include "bfgs.h"
#include "callback_handler.h"
#include "deriv_check.h"
#include "dual_estimation.h"
#include "linesearch.h"
#include "gauss_newton.h"
#include "merit.h"
#include "newton.h"
#include "parametric.h"
#include "polish.h"
#include "problem_scaling.h"
#include "soc.h"
#include "sr1.h"

#include "cauchy/cauchy.h"
#include "lp/lpi.h"
#include "step/step_rule.h"
#include "sparse/sparse_factorization.h"

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

    SleqpStepRule* step_rule;

    SleqpIterate* iterate;

    SleqpIterate* trial_iterate;

    SleqpSparseVec* original_violation;

    SleqpLPi* lp_interface;

    SleqpCauchy* cauchy_data;

    SleqpSparseVec* cauchy_direction;

    SleqpSparseVec* cauchy_step;

    SleqpSparseVec* cauchy_hessian_step;

    double cauchy_step_length;

    SleqpSparseVec* multipliers;

    SleqpWorkingStep* working_step;

    SleqpNewtonData* newton_data;

    SleqpGaussNewtonSolver* gauss_newton_solver;

    SleqpSparseVec* newton_step;

    SleqpSparseVec* newton_hessian_step;

    SleqpSparseVec* trial_step;

    SLEQP_STEPTYPE last_step_type;

    SleqpSparseVec* initial_trial_point;

    SleqpSparseFactorization* factorization;

    SleqpAugJac* aug_jac;

    SleqpDualEstimation* estimation_data;
    SleqpSparseVec* estimation_residuals;

    SleqpMeritData* merit_data;

    SleqpLineSearchData* linesearch;

    SleqpPolishing* polishing;

    SleqpParametricSolver* parametric_solver;
    SleqpWorkingSet* parametric_original_working_set;

    SleqpCallbackHandler** callback_handlers;

    // Primal / dual step lengths

    SleqpSparseVec* primal_diff;

    double primal_diff_norm;

    SleqpSparseVec* cons_dual_diff;

    SleqpSparseVec* vars_dual_diff;

    double dual_diff_norm;

    double current_merit_value;

    // SOC related
    SleqpSOC* soc_data;

    SleqpSparseVec* soc_step;

    double* dense_cache;

    // residuum

    double slackness_residuum;

    double stationarity_residuum;

    double feasibility_residuum;

    // BFGS related

    SleqpBFGS* bfgs_data;

    // SR1 related

    SleqpSR1* sr1_data;

    // parameters, adjusted throughout...

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
                                                    double cauchy_step_length,
                                                    double eps,
                                                    double* lp_trust_radius);

  SLEQP_RETCODE sleqp_solver_compute_trial_point_simple(SleqpSolver* solver,
                                                        double* cauchy_merit_value,
                                                        bool quadratic_model,
                                                        bool* full_step);

  SLEQP_RETCODE sleqp_solver_compute_trial_point_newton(SleqpSolver* solver,
                                                        double* trial_merit_value,
                                                        bool* full_step);

  SLEQP_RETCODE sleqp_solver_compute_trial_point(SleqpSolver* solver,
                                                 double* trial_merit_value,
                                                 bool* full_step,
                                                 bool* reject);

  SLEQP_RETCODE sleqp_solver_compute_trial_point_det(SleqpSolver* solver,
                                                     double* trial_merit_value,
                                                     bool* full_step);

  SLEQP_RETCODE sleqp_solver_compute_trial_point_soc(SleqpSolver* solver,
                                                     bool* reject);

  SLEQP_RETCODE sleqp_solver_compute_cauchy_step(SleqpSolver* solver,
                                                 double* cauchy_merit_value,
                                                 bool quadratic_model,
                                                 bool* full_step);

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
