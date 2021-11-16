#ifndef SLEQP_TRIAL_POINT_H
#define SLEQP_TRIAL_POINT_H

#include "dual_estimation.h"
#include "eqp.h"
#include "iterate.h"
#include "linesearch.h"
#include "merit.h"
#include "options.h"
#include "parametric.h"
#include "problem.h"
#include "soc.h"
#include "working_step.h"

#include "cauchy/cauchy.h"
#include "factorization/factorization.h"
#include "lp/lpi.h"

typedef struct
{
  int refcount;
  SleqpProblem* problem;

  SleqpParams* params;
  SleqpOptions* options;

  SleqpSparseVec* cauchy_direction;

  SleqpSparseVec* cauchy_step;
  SleqpSparseVec* cauchy_hessian_step;

  SleqpSparseVec* estimation_residuals;

  SleqpSparseVec* newton_step;
  SleqpSparseVec* newton_hessian_step;

  SleqpSparseVec* soc_step;

  SleqpSparseVec* trial_step;

  SleqpSparseVec* multipliers;

  SleqpSparseVec* initial_trial_point;

  SleqpMerit* merit;

  SleqpIterate* iterate;

  SleqpLPi* lp_interface;

  SleqpCauchy* cauchy_data;

  SleqpDualEstimation* estimation_data;

  SleqpFactorization* factorization;

  SleqpAugJac* aug_jac;

  SleqpLineSearch* linesearch;

  SleqpWorkingStep* working_step;

  SleqpEQPSolver* eqp_solver;

  SleqpSOC* soc_data;

  SleqpParametricSolver* parametric_solver;
  SleqpWorkingSet* parametric_original_working_set;

  double* dense_cache;

  SleqpTimer* elapsed_timer;

  double penalty_parameter;

  double lp_trust_radius;
  double trust_radius;

  bool locally_infeasible;

  double current_merit_value;

  double time_limit;

} SleqpTrialPointSolver;

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_trial_point_solver_create(SleqpTrialPointSolver** star,
                                SleqpProblem* problem,
                                SleqpParams* params,
                                SleqpOptions* options);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_trial_point_solver_set_iterate(SleqpTrialPointSolver* solver,
                                     SleqpIterate* iterate);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_trial_point_solver_set_time_limit(SleqpTrialPointSolver* solver,
                                        double time_limit);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_trial_point_solver_set_trust_radius(SleqpTrialPointSolver* solver,
                                          double trust_radius);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_trial_point_solver_set_lp_trust_radius(SleqpTrialPointSolver* solver,
                                             double lp_trust_radius);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_trial_point_solver_set_penalty(SleqpTrialPointSolver* solver,
                                     double penalty_parameter);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_trial_point_solver_get_penalty(SleqpTrialPointSolver* solver,
                                     double* penalty_parameter);

bool
sleqp_trial_point_solver_locally_infeasible(SleqpTrialPointSolver* solver);

SleqpSparseVec*
sleqp_trial_point_solver_get_multipliers(SleqpTrialPointSolver* solver);

SleqpSparseVec*
sleqp_trial_point_solver_get_cauchy_step(SleqpTrialPointSolver* solver);

SleqpSparseVec*
sleqp_trial_point_solver_get_trial_step(SleqpTrialPointSolver* solver);

SleqpSparseVec*
sleqp_trial_point_solver_get_soc_step(SleqpTrialPointSolver* solver);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_trial_point_solver_get_rayleigh(SleqpTrialPointSolver* solver,
                                      double* min_rayleigh,
                                      double* max_rayleigh);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_trial_point_solver_print_stats(SleqpTrialPointSolver* solver,
                                     double elapsed_seconds);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_trial_point_solver_compute_cauchy_step(SleqpTrialPointSolver* solver,
                                             double* cauchy_merit_value,
                                             bool quadratic_model,
                                             bool* full_step);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_trial_point_solver_compute_trial_point(SleqpTrialPointSolver* solver,
                                             SleqpIterate* trial_iterate,
                                             double* trial_merit_value,
                                             bool* full_step,
                                             bool* reject);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_trial_point_solver_compute_trial_point_soc(SleqpTrialPointSolver* solver,
                                                 SleqpIterate* trial_iterate,
                                                 bool* reject);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_trial_point_solver_capture(SleqpTrialPointSolver* solver);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_trial_point_solver_release(SleqpTrialPointSolver** star);

#endif /* SLEQP_TRIAL_POINT_H */
