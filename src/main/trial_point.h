#ifndef SLEQP_TRIAL_POINT_H
#define SLEQP_TRIAL_POINT_H

#include "direction.h"
#include "dual_estimation/dual_estimation.h"
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
#include "fact/fact.h"
#include "lp/lpi.h"

typedef struct
{
  int refcount;
  SleqpProblem* problem;

  SleqpParams* params;
  SleqpOptions* options;

  SleqpVec* lp_step;

  SleqpDirection* cauchy_direction;
  SleqpVec* estimation_residuals;

  SleqpDirection* newton_direction;

  SleqpDirection* soc_direction;

  SleqpDirection* trial_direction;

  SleqpVec* multipliers;

  SleqpVec* initial_trial_point;

  SleqpMerit* merit;

  SleqpIterate* iterate;

  SleqpCauchy* cauchy_data;

  SleqpDualEstimation* estimation_data;

  SleqpFact* fact;

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

  double feasibility_residuum;
  bool allow_global_reset;
  bool performed_global_reset;

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
sleqp_trial_point_solver_set_penalty_info(SleqpTrialPointSolver* solver,
                                          double feas_res,
                                          bool allow_global_reset);

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
sleqp_trial_point_solver_penalty(SleqpTrialPointSolver* solver,
                                 double* penalty_parameter);

bool
sleqp_trial_point_solver_locally_infeasible(SleqpTrialPointSolver* solver);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_trial_point_solver_penalty_info(SleqpTrialPointSolver* solver,
                                      bool* performed_global_reset);

SleqpVec*
sleqp_trial_point_solver_multipliers(SleqpTrialPointSolver* solver);

SleqpVec*
sleqp_trial_point_solver_cauchy_step(SleqpTrialPointSolver* solver);

SleqpVec*
sleqp_trial_point_solver_trial_step(SleqpTrialPointSolver* solver);

SleqpVec*
sleqp_trial_point_solver_soc_step(SleqpTrialPointSolver* solver);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_trial_point_solver_rayleigh(SleqpTrialPointSolver* solver,
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
                                             bool* failed_eqp_step,
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
