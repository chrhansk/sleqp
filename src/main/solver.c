#include "solver.h"

#include <assert.h>
#include <math.h>
#include <threads.h>

#include "aug_jacobian.h"
#include "bfgs.h"
#include "callback_handler.h"
#include "cauchy.h"
#include "cmp.h"
#include "defs.h"
#include "deriv_check.h"
#include "dual_estimation.h"
#include "fail.h"
#include "feas.h"
#include "iterate.h"
#include "linesearch.h"
#include "mem.h"
#include "merit.h"
#include "newton.h"
#include "parametric.h"
#include "penalty.h"
#include "problem_scaling.h"
#include "scale.h"
#include "soc.h"
#include "sr1.h"
#include "timer.h"
#include "util.h"

#include "lp/lpi.h"
#include "preprocessor/preprocessor.h"

#define INFO_BUF_SIZE 100
#define SOLVER_INFO_BUF_SIZE 400

thread_local char lps_info[INFO_BUF_SIZE];
thread_local char fact_info[INFO_BUF_SIZE];
thread_local char solver_info[SOLVER_INFO_BUF_SIZE];

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

  SleqpNewtonData* newton_data;

  SleqpSparseVec* newton_step;

  SleqpSparseVec* newton_hessian_step;

  SleqpSparseVec* trial_step;

  SLEQP_STEPTYPE last_step_type;

  SleqpSparseVec* initial_trial_point;

  SleqpSparseFactorization* factorization;

  SleqpAugJacobian* aug_jacobian;

  SleqpDualEstimation* estimation_data;
  SleqpSparseVec* estimation_residuals;

  SleqpMeritData* merit_data;

  SleqpLineSearchData* linesearch;

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

  SleqpSparseVec* soc_trial_point;

  double soc_step_norm;

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

static double remaining_time(SleqpSolver* solver)
{
  double time_limit = solver->time_limit;

  if(time_limit != SLEQP_NONE)
  {
    double remaining_time = time_limit - sleqp_timer_get_ttl(solver->elapsed_timer);

    remaining_time = SLEQP_MAX(remaining_time, 0.);

    return remaining_time;
  }

  return SLEQP_NONE;
}

static SLEQP_RETCODE set_residuum(SleqpSolver* solver)
{
  SLEQP_CALL(sleqp_iterate_slackness_residuum(solver->problem,
                                              solver->iterate,
                                              &solver->slackness_residuum));

  SLEQP_CALL(sleqp_iterate_feasibility_residuum(solver->problem,
                                                solver->iterate,
                                                &solver->feasibility_residuum));

  SLEQP_CALL(sleqp_iterate_stationarity_residuum(solver->problem,
                                                 solver->iterate,
                                                 solver->dense_cache,
                                                 &solver->stationarity_residuum));

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE solver_convert_primal(SleqpSolver* solver,
                                    const SleqpSparseVec* source,
                                    SleqpSparseVec* target)
{
  assert(source->dim == sleqp_problem_num_variables(solver->original_problem));
  assert(target->dim == sleqp_problem_num_variables(solver->problem));

  SLEQP_CALL(sleqp_sparse_vector_copy(source, solver->scaled_primal));

  if(solver->scaling_data)
  {
    SLEQP_CALL(sleqp_scale_point(solver->scaling_data,
                                 solver->scaled_primal));
  }

  if(solver->preprocessor)
  {
    SLEQP_CALL(sleqp_preprocessor_transform_primal(solver->preprocessor,
                                                   solver->scaled_primal,
                                                   target));
  }
  else
  {
    SLEQP_CALL(sleqp_sparse_vector_copy(solver->scaled_primal, target));
  }


  return SLEQP_OKAY;
}

static
SLEQP_RETCODE solver_restore_iterate(SleqpSolver* solver,
                                     const SleqpIterate* source,
                                     SleqpIterate* target)
{
  if(solver->preprocessor)
  {
    SLEQP_CALL(sleqp_preprocessor_restore_iterate(solver->preprocessor,
                                                  source,
                                                  solver->scaled_iterate));
  }
  else
  {
    SLEQP_CALL(sleqp_iterate_copy(source, solver->scaled_iterate));
  }

  SLEQP_CALL(sleqp_iterate_copy(solver->scaled_iterate, target));

  if(solver->scaling_data)
  {
    SLEQP_CALL(sleqp_unscale_iterate(solver->scaling_data, target));
  }

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE restore_original_iterate(SleqpSolver* solver)
{
  if(solver->restore_original_iterate)
  {
    SLEQP_CALL(solver_restore_iterate(solver,
                                      solver->iterate,
                                      solver->original_iterate));

    solver->restore_original_iterate = false;
  }

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE solver_create_problem(SleqpSolver* solver,
                                    SleqpProblem* problem)
{
  SleqpProblem* scaled_problem;

  SleqpParams* params = solver->params;
  SleqpOptions* options = solver->options;

  solver->original_problem = problem;
  SLEQP_CALL(sleqp_problem_capture(solver->original_problem));

  if(solver->scaling_data)
  {
    SLEQP_CALL(sleqp_problem_scaling_create(&solver->problem_scaling,
                                            solver->scaling_data,
                                            problem,
                                            params,
                                            options));

    SLEQP_CALL(sleqp_problem_scaling_flush(solver->problem_scaling));

    scaled_problem = sleqp_problem_scaling_get_problem(solver->problem_scaling);
  }
  else
  {
    scaled_problem = problem;
  }

  SleqpFunc* func = sleqp_problem_func(scaled_problem);

  {
    const SLEQP_HESSIAN_EVAL hessian_eval = sleqp_options_get_int(options,
                                                                  SLEQP_OPTION_INT_HESSIAN_EVAL);

    if(hessian_eval == SLEQP_HESSIAN_EVAL_SIMPLE_BFGS ||
       hessian_eval == SLEQP_HESSIAN_EVAL_DAMPED_BFGS)
    {
      SLEQP_CALL(sleqp_bfgs_create(&solver->bfgs_data,
                                   func,
                                   params,
                                   options));

      func = sleqp_bfgs_get_func(solver->bfgs_data);
    }

    if(hessian_eval == SLEQP_HESSIAN_EVAL_SR1)
    {
      SLEQP_CALL(sleqp_sr1_create(&solver->sr1_data,
                                  func,
                                  params,
                                  options));

      func = sleqp_sr1_get_func(solver->sr1_data);
    }

    SLEQP_CALL(sleqp_problem_create(&solver->scaled_problem,
                                    func,
                                    params,
                                    sleqp_problem_var_lb(scaled_problem),
                                    sleqp_problem_var_ub(scaled_problem),
                                    sleqp_problem_general_lb(scaled_problem),
                                    sleqp_problem_general_ub(scaled_problem),
                                    sleqp_problem_linear_coeffs(scaled_problem),
                                    sleqp_problem_linear_lb(scaled_problem),
                                    sleqp_problem_linear_ub(scaled_problem)));
  }

  const bool enable_preprocesor = sleqp_options_get_bool(solver->options,
                                                         SLEQP_OPTION_BOOL_ENABLE_PREPROCESSOR);

  if(enable_preprocesor)
  {
    SLEQP_CALL(sleqp_preprocessor_create(&solver->preprocessor,
                                         solver->scaled_problem,
                                         solver->params));

    const SLEQP_PREPROCESSING_RESULT preprocessing_result = sleqp_preprocessor_result(solver->preprocessor);

    if(preprocessing_result == SLEQP_PREPROCESSING_RESULT_FAILURE)
    {
      solver->problem = solver->scaled_problem;
    }
    else
    {
      solver->problem = sleqp_preprocessor_transformed_problem(solver->preprocessor);

      if(preprocessing_result == SLEQP_PREPROCESSING_RESULT_INFEASIBLE)
      {
        // ...
      }
    }
  }
  else
  {
    solver->problem = solver->scaled_problem;
  }

  SLEQP_CALL(sleqp_problem_capture(solver->problem));


  return SLEQP_OKAY;
}

static
SLEQP_RETCODE solver_create_iterates(SleqpSolver* solver,
                                     SleqpSparseVec* primal)
{
  const double zero_eps = sleqp_params_get(solver->params,
                                           SLEQP_PARAM_ZERO_EPS);

  SLEQP_CALL(solver_convert_primal(solver,
                                   primal,
                                   solver->primal));

  SLEQP_CALL(sleqp_iterate_create(&solver->iterate,
                                  solver->problem,
                                  solver->primal));

  SLEQP_CALL(sleqp_sparse_vector_clip(solver->primal,
                                      sleqp_problem_var_lb(solver->problem),
                                      sleqp_problem_var_ub(solver->problem),
                                      zero_eps,
                                      sleqp_iterate_get_primal(solver->iterate)));

  SLEQP_CALL(sleqp_iterate_create(&solver->trial_iterate,
                                  solver->problem,
                                  sleqp_iterate_get_primal(solver->iterate)));

  SLEQP_CALL(sleqp_iterate_create(&solver->scaled_iterate,
                                  solver->scaled_problem,
                                  primal));

  if(solver->scaling_data || solver->preprocessor)
  {
    SLEQP_CALL(sleqp_iterate_create(&solver->original_iterate,
                                    solver->original_problem,
                                    primal));

    solver->restore_original_iterate = true;
  }
  else
  {
    solver->original_iterate = solver->iterate;
    solver->restore_original_iterate = false;
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solver_create(SleqpSolver** star,
                                  SleqpProblem* problem,
                                  SleqpParams* params,
                                  SleqpOptions* options,
                                  SleqpSparseVec* primal,
                                  SleqpScaling* scaling_data)
{
  assert(sleqp_sparse_vector_is_valid(primal));

  SLEQP_CALL(sleqp_malloc(star));

  SleqpSolver* solver = *star;

  *solver = (SleqpSolver){0};

  solver->refcount = 1;

  const int num_original_variables = sleqp_problem_num_variables(problem);

  SLEQP_CALL(sleqp_timer_create(&solver->elapsed_timer));

  if(scaling_data)
  {
    SLEQP_CALL(sleqp_scaling_capture(scaling_data));
    solver->scaling_data = scaling_data;
  }

  SLEQP_CALL(sleqp_params_capture(params));
  solver->params = params;

  SLEQP_CALL(sleqp_options_capture(options));
  solver->options = options;

  SLEQP_CALL(solver_create_problem(solver,
                                   problem));

  const int num_variables = sleqp_problem_num_variables(solver->problem);
  const int num_constraints = sleqp_problem_num_constraints(solver->problem);

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->scaled_primal,
                                              num_original_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->primal,
                                              num_variables));

  SLEQP_CALL(solver_create_iterates(solver,
                                    primal));

  solver->iteration = 0;
  solver->elapsed_seconds = 0.;

  SLEQP_CALL(sleqp_deriv_checker_create(&solver->deriv_check,
                                        solver->problem,
                                        params));

  const int num_lp_variables = num_variables + 2*num_constraints;
  const int num_lp_constraints = num_constraints;

  SLEQP_CALL(sleqp_lpi_create_default_interface(&solver->lp_interface,
                                                num_lp_variables,
                                                num_lp_constraints,
                                                params,
                                                options));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->original_violation,
                                              num_constraints));

  SLEQP_CALL(sleqp_cauchy_create(&solver->cauchy_data,
                                 solver->problem,
                                 params,
                                 options,
                                 solver->lp_interface));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->cauchy_direction,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->cauchy_step,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->cauchy_hessian_step,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->multipliers,
                                              num_constraints));

  SLEQP_CALL(sleqp_newton_data_create(&solver->newton_data,
                                      solver->problem,
                                      params,
                                      options));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->newton_step,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->newton_hessian_step,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->trial_step,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->initial_trial_point,
                                              num_variables));

  // create sparse factorization

  SLEQP_CALL(sleqp_sparse_factorization_create_default(&solver->factorization,
                                                       params));

  SLEQP_CALL(sleqp_aug_jacobian_create(&solver->aug_jacobian,
                                       solver->problem,
                                       params,
                                       solver->factorization));

  SLEQP_CALL(sleqp_dual_estimation_create(&solver->estimation_data,
                                          solver->problem));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->estimation_residuals,
                                              num_variables));

  SLEQP_CALL(sleqp_merit_data_create(&solver->merit_data,
                                     solver->problem,
                                     params));

  SLEQP_CALL(sleqp_linesearch_create(&solver->linesearch,
                                     solver->problem,
                                     params,
                                     solver->merit_data));

  SLEQP_PARAMETRIC_CAUCHY parametric_cauchy = sleqp_options_get_int(solver->options,
                                                                    SLEQP_OPTION_INT_PARAMETRIC_CAUCHY);

  if(parametric_cauchy != SLEQP_PARAMETRIC_CAUCHY_DISABLED)
  {
    SLEQP_CALL(sleqp_parametric_solver_create(&solver->parametric_solver,
                                              solver->problem,
                                              solver->params,
                                              solver->options,
                                              solver->merit_data,
                                              solver->linesearch));

    SLEQP_CALL(sleqp_working_set_create(&solver->parametric_original_working_set, solver->problem));
  }

  SLEQP_CALL(sleqp_alloc_array(&solver->callback_handlers,
                               SLEQP_SOLVER_NUM_EVENTS));

  for(int i = 0; i < SLEQP_SOLVER_NUM_EVENTS; ++i)
  {
    SLEQP_CALL(sleqp_callback_handler_create(solver->callback_handlers + i));
  }

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->primal_diff,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->cons_dual_diff,
                                              num_constraints));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->vars_dual_diff,
                                              num_variables));

  SLEQP_CALL(sleqp_soc_data_create(&solver->soc_data,
                                   solver->problem,
                                   params));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->soc_trial_point,
                                              num_variables));

  SLEQP_CALL(sleqp_alloc_array(&solver->dense_cache, SLEQP_MAX(num_variables, num_constraints)));

  SLEQP_CALL(sleqp_solver_reset(solver));

  solver->time_limit = SLEQP_NONE;

  solver->abort_next = false;

  sleqp_log_debug("%s", sleqp_solver_info(solver));

  return SLEQP_OKAY;
}

const char* sleqp_solver_info(const SleqpSolver* solver)
{
  snprintf(lps_info,
           INFO_BUF_SIZE,
           "%s %s",
           sleqp_lpi_get_name(solver->lp_interface),
           sleqp_lpi_get_version(solver->lp_interface));

  snprintf(fact_info,
           INFO_BUF_SIZE,
           "%s %s",
           sleqp_sparse_factorization_get_name(solver->factorization),
           sleqp_sparse_factorization_get_version(solver->factorization));

  snprintf(solver_info,
           SOLVER_INFO_BUF_SIZE,
           "Sleqp version %s [LP solver: %s] [factorization: %s] [GitHash %s]",
           SLEQP_VERSION,
           lps_info,
           fact_info,
           SLEQP_GIT_COMMIT_HASH);

  return solver_info;
}

static SLEQP_RETCODE update_lp_trust_radius(bool trial_step_accepted,
                                            double trial_step_infnorm,
                                            double cauchy_step_infnorm,
                                            double cauchy_step_length,
                                            double eps,
                                            double* lp_trust_radius)
{
  if(trial_step_accepted)
  {
    double norm_increase_factor = 1.2;

    trial_step_infnorm *= norm_increase_factor;
    cauchy_step_infnorm *= norm_increase_factor;

    double scaled_trust_radius = .1 * (*lp_trust_radius);

    double update_lhs = SLEQP_MAX(trial_step_infnorm,
                                  cauchy_step_infnorm);

    update_lhs = SLEQP_MAX(update_lhs, scaled_trust_radius);

    if(sleqp_is_eq(cauchy_step_length, 1., eps))
    {
      (*lp_trust_radius) *= 7.;
    }

    *lp_trust_radius = SLEQP_MIN(update_lhs, *lp_trust_radius);

  }
  else
  {
    double half_norm = .5 * trial_step_infnorm;
    double small_radius = .1 * (*lp_trust_radius);

    double reduced_radius = SLEQP_MAX(half_norm, small_radius);

    *lp_trust_radius = SLEQP_MIN(reduced_radius, *lp_trust_radius);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE update_trust_radius(SleqpSolver* solver,
                                         double reduction_ratio,
                                         bool trial_step_accepted,
                                         double direction_norm)
{
  const double eps = sleqp_params_get(solver->params, SLEQP_PARAM_EPS);

  double* trust_radius = &(solver->trust_radius);

  if(reduction_ratio >= 0.9)
  {
    *trust_radius = SLEQP_MAX(*trust_radius,
                              7*direction_norm);
  }
  else if(reduction_ratio >= 0.3)
  {
    *trust_radius = SLEQP_MAX(*trust_radius,
                              2*direction_norm);
  }
  else if(trial_step_accepted)
  {
    // stays the same
  }
  else
  {
    // filter out very small steps
    if(sleqp_is_zero(direction_norm, eps))
    {
      *trust_radius *= .5;
    }
    else
    {
      *trust_radius = SLEQP_MIN(.5 * (*trust_radius),
                                .5 * direction_norm);
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE estimate_dual_values(SleqpSolver* solver,
                                          SleqpIterate* iterate)
{
  SleqpOptions* options = solver->options;

  SLEQP_DUAL_ESTIMATION_TYPE estimation_type = sleqp_options_get_int(options,
                                                                     SLEQP_OPTION_INT_DUAL_ESTIMATION_TYPE);

  if(estimation_type == SLEQP_DUAL_ESTIMATION_TYPE_LSQ)
  {
    SLEQP_CALL(sleqp_dual_estimation_compute(solver->estimation_data,
                                             iterate,
                                             solver->estimation_residuals,
                                             solver->aug_jacobian));
  }
  else
  {
    assert(estimation_type == SLEQP_DUAL_ESTIMATION_TYPE_LP);

    SLEQP_CALL(sleqp_cauchy_get_dual_estimation(solver->cauchy_data,
                                                iterate));
  }

  SLEQP_CALL(sleqp_newton_add_violated_multipliers(solver->newton_data,
                                                   solver->multipliers));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE compute_cauchy_step_parametric(SleqpSolver* solver,
                                                    double* cauchy_merit_value,
                                                    bool* full_step)
{
  SleqpIterate* iterate = solver->iterate;

  const double eps = sleqp_params_get(solver->params, SLEQP_PARAM_EPS);

  SLEQP_NUM_ASSERT_PARAM(eps);

  {
    SLEQP_CALL(sleqp_cauchy_set_iterate(solver->cauchy_data,
                                        iterate,
                                        solver->lp_trust_radius));

    SLEQP_CALL(sleqp_cauchy_solve(solver->cauchy_data,
                                  sleqp_iterate_get_func_grad(iterate),
                                  solver->penalty_parameter,
                                  SLEQP_CAUCHY_OBJECTIVE_TYPE_DEFAULT));

    double objective_value;

    SLEQP_CALL(sleqp_cauchy_get_objective_value(solver->cauchy_data,
                                                &objective_value));

    const double reduction = solver->current_merit_value - objective_value;

    sleqp_assert_is_geq(reduction, 0., eps);

    const double criticality_bound = reduction / SLEQP_MIN(solver->lp_trust_radius, 1.);

    // Bound on the criticality measure used in
    // "On the Convergence of Successive Linear Programming Algorithms"
    sleqp_log_debug("Criticality bound: %g", criticality_bound);


    SLEQP_CALL(sleqp_cauchy_get_working_set(solver->cauchy_data,
                                            iterate));

    SLEQP_CALL(sleqp_aug_jacobian_set_iterate(solver->aug_jacobian,
                                              iterate));

    SLEQP_CALL(sleqp_newton_set_iterate(solver->newton_data,
                                        iterate,
                                        solver->aug_jacobian,
                                        solver->trust_radius,
                                        solver->penalty_parameter));

    SLEQP_CALL(estimate_dual_values(solver, iterate));

    SLEQP_CALL(sleqp_cauchy_get_direction(solver->cauchy_data,
                                          solver->cauchy_direction));
  }

  {
    SLEQP_CALL(sleqp_parametric_solver_set_penalty(solver->parametric_solver,
                                                   solver->penalty_parameter));

    SLEQP_CALL(sleqp_parametric_solver_solve(solver->parametric_solver,
                                             iterate,
                                             solver->cauchy_data,
                                             solver->cauchy_step,
                                             solver->cauchy_hessian_step,
                                             solver->multipliers,
                                             &(solver->lp_trust_radius),
                                             cauchy_merit_value));
  }

  SLEQP_CALL(sleqp_working_set_copy(sleqp_iterate_get_working_set(iterate),
                                    solver->parametric_original_working_set));

  SLEQP_CALL(sleqp_cauchy_get_working_set(solver->cauchy_data,
                                          iterate));

  // Reconstruct the augmented Jacobian if required
  if(!sleqp_working_set_eq(solver->parametric_original_working_set,
                           sleqp_iterate_get_working_set(iterate)))
  {
    SLEQP_CALL(sleqp_aug_jacobian_set_iterate(solver->aug_jacobian,
                                              iterate));
  }

  SLEQP_CALL(sleqp_linesearch_set_iterate(solver->linesearch,
                                          iterate,
                                          solver->penalty_parameter,
                                          solver->trust_radius));

  SLEQP_CALL(sleqp_newton_set_iterate(solver->newton_data,
                                      iterate,
                                      solver->aug_jacobian,
                                      solver->trust_radius,
                                      solver->penalty_parameter));

#if !defined(NDEBUG)

  {
    double actual_quadratic_merit_value;

    double func_dual = 1.;

    SLEQP_CALL(sleqp_merit_quadratic(solver->merit_data,
                                     iterate,
                                     &func_dual,
                                     solver->cauchy_step,
                                     solver->multipliers,
                                     solver->penalty_parameter,
                                     &actual_quadratic_merit_value));

    sleqp_assert_is_eq(*cauchy_merit_value,
                       actual_quadratic_merit_value,
                       eps);
  }

#endif

  (*full_step) = true;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE compute_cauchy_step(SleqpSolver* solver,
                                         double* cauchy_merit_value,
                                         bool quadratic_model,
                                         bool* full_step)
{
  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;

  const double eps = sleqp_params_get(solver->params,
                                      SLEQP_PARAM_EPS);

  SLEQP_NUM_ASSERT_PARAM(eps);

  const double zero_eps = sleqp_params_get(solver->params,
                                           SLEQP_PARAM_ZERO_EPS);

  const double one = 1.;

  // compute Cauchy direction / step and dual estimation
  {
    SLEQP_CALL(sleqp_cauchy_set_iterate(solver->cauchy_data,
                                        iterate,
                                        solver->lp_trust_radius));

    SLEQP_CALL(sleqp_cauchy_solve(solver->cauchy_data,
                                  sleqp_iterate_get_func_grad(iterate),
                                  solver->penalty_parameter,
                                  SLEQP_CAUCHY_OBJECTIVE_TYPE_DEFAULT));

    double objective_value;

    SLEQP_CALL(sleqp_cauchy_get_objective_value(solver->cauchy_data,
                                                &objective_value));

    const double reduction = solver->current_merit_value - objective_value;

    sleqp_assert_is_geq(reduction, 0., eps);

    const double criticality_bound = reduction / SLEQP_MIN(solver->lp_trust_radius, 1.);

    // Bound on the criticality measure used in
    // "On the Convergence of Successive Linear Programming Algorithms"
    sleqp_log_debug("Criticality bound: %g", criticality_bound);


    SLEQP_CALL(sleqp_cauchy_get_working_set(solver->cauchy_data,
                                            iterate));

    SLEQP_CALL(sleqp_aug_jacobian_set_iterate(solver->aug_jacobian,
                                              iterate));

    SLEQP_CALL(sleqp_newton_set_iterate(solver->newton_data,
                                        iterate,
                                        solver->aug_jacobian,
                                        solver->trust_radius,
                                        solver->penalty_parameter));

    SLEQP_CALL(estimate_dual_values(solver, iterate));

    SLEQP_CALL(sleqp_cauchy_get_direction(solver->cauchy_data,
                                          solver->cauchy_direction));

#if !defined(NDEBUG)

    {
      bool in_working_set = false;

      SLEQP_CALL(sleqp_direction_in_working_set(problem,
                                                iterate,
                                                solver->cauchy_direction,
                                                solver->dense_cache,
                                                eps,
                                                &in_working_set));

      assert(in_working_set);
    }

#endif

    SLEQP_CALL(sleqp_sparse_vector_copy(solver->cauchy_direction,
                                        solver->cauchy_step));

    if(!quadratic_model)
    {
      (*full_step) = true;

      solver->cauchy_step_length = 1.;

      return SLEQP_OKAY;
    }

    SLEQP_CALL(sleqp_problem_hess_prod(problem,
                                       &one,
                                       solver->cauchy_direction,
                                       solver->multipliers,
                                       solver->cauchy_hessian_step));

    SLEQP_CALL(sleqp_linesearch_set_iterate(solver->linesearch,
                                            iterate,
                                            solver->penalty_parameter,
                                            solver->trust_radius));

    SLEQP_CALL(sleqp_linesearch_cauchy_step(solver->linesearch,
                                            solver->cauchy_step,
                                            solver->multipliers,
                                            solver->cauchy_hessian_step,
                                            &solver->cauchy_step_length,
                                            cauchy_merit_value));

#if !defined(NDEBUG)

    {
      double actual_quadratic_merit_value;

      double func_dual = 1.;

      SLEQP_CALL(sleqp_merit_quadratic(solver->merit_data,
                                       iterate,
                                       &func_dual,
                                       solver->cauchy_step,
                                       solver->multipliers,
                                       solver->penalty_parameter,
                                       &actual_quadratic_merit_value));

      sleqp_assert_is_eq(*cauchy_merit_value,
                         actual_quadratic_merit_value,
                         eps);
    }

#endif

    (*full_step) = sleqp_is_eq(solver->cauchy_step_length,
                               1.,
                               zero_eps);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE compute_trial_point_simple(SleqpSolver* solver,
                                                double* cauchy_merit_value,
                                                bool quadratic_model,
                                                bool* full_step)
{
  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;

  const double eps = sleqp_params_get(solver->params,
                                      SLEQP_PARAM_EPS);

  SLEQP_NUM_ASSERT_PARAM(eps);

  const double zero_eps = sleqp_params_get(solver->params,
                                           SLEQP_PARAM_ZERO_EPS);

  SLEQP_PARAMETRIC_CAUCHY parametric_cauchy = sleqp_options_get_int(solver->options,
                                                                    SLEQP_OPTION_INT_PARAMETRIC_CAUCHY);

  if(parametric_cauchy != SLEQP_PARAMETRIC_CAUCHY_DISABLED)
  {
    SLEQP_CALL(compute_cauchy_step_parametric(solver,
                                              cauchy_merit_value,
                                              full_step));
  }
  else
  {
    SLEQP_CALL(compute_cauchy_step(solver,
                                   cauchy_merit_value,
                                   quadratic_model,
                                   full_step));
  }

  const SleqpSparseVec* trial_step = solver->cauchy_step;

  // Compute quadratic merit value
  {
    SLEQP_CALL(sleqp_merit_linear(solver->merit_data,
                                  solver->iterate,
                                  trial_step,
                                  solver->penalty_parameter,
                                  cauchy_merit_value));
  }

  if(quadratic_model)
  {
    double hessian_prod;

    SLEQP_CALL(sleqp_sparse_vector_dot(trial_step,
                                       solver->cauchy_hessian_step,
                                       &hessian_prod));

    (*cauchy_merit_value) += .5 * hessian_prod;

#if !defined(NDEBUG)

    {
      double actual_quadratic_merit_value;

      double func_dual = 1.;

      SLEQP_CALL(sleqp_merit_quadratic(solver->merit_data,
                                       iterate,
                                       &func_dual,
                                       solver->cauchy_step,
                                       solver->multipliers,
                                       solver->penalty_parameter,
                                       &actual_quadratic_merit_value));

      sleqp_assert_is_eq(*cauchy_merit_value,
                         actual_quadratic_merit_value,
                         eps);
    }

#endif

  }

  SLEQP_CALL(sleqp_sparse_vector_add(sleqp_iterate_get_primal(iterate),
                                     solver->cauchy_step,
                                     zero_eps,
                                     solver->initial_trial_point));

  SLEQP_CALL(sleqp_sparse_vector_clip(solver->initial_trial_point,
                                      sleqp_problem_var_lb(problem),
                                      sleqp_problem_var_ub(problem),
                                      zero_eps,
                                      sleqp_iterate_get_primal(solver->trial_iterate)));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE compute_trial_point_newton(SleqpSolver* solver,
                                                double* trial_merit_value,
                                                bool* full_step)
{
  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;

  const double eps = sleqp_params_get(solver->params, SLEQP_PARAM_EPS);

  SLEQP_NUM_ASSERT_PARAM(eps);

  const double zero_eps = sleqp_params_get(solver->params,
                                           SLEQP_PARAM_ZERO_EPS);

  const double one = 1.;

  double cauchy_merit_value;

  SLEQP_PARAMETRIC_CAUCHY parametric_cauchy = sleqp_options_get_int(solver->options,
                                                                    SLEQP_OPTION_INT_PARAMETRIC_CAUCHY);

  if(parametric_cauchy != SLEQP_PARAMETRIC_CAUCHY_DISABLED)
  {
    SLEQP_CALL(compute_cauchy_step_parametric(solver,
                                              &cauchy_merit_value,
                                              full_step));
  }
  else
  {
    SLEQP_CALL(compute_cauchy_step(solver,
                                   &cauchy_merit_value,
                                   true,
                                   full_step));
  }

  // compute Newton step
  {
    SLEQP_CALL(sleqp_newton_set_time_limit(solver->newton_data,
                                           remaining_time(solver)));

    SLEQP_CALL(sleqp_newton_compute_step(solver->newton_data,
                                         solver->multipliers,
                                         solver->newton_step));

    SLEQP_CALL(sleqp_problem_hess_prod(problem,
                                       &one,
                                       solver->newton_step,
                                       solver->multipliers,
                                       solver->newton_hessian_step));
  }

  {
    SLEQP_LINESEARCH lineserach = sleqp_options_get_int(solver->options,
                                                        SLEQP_OPTION_INT_LINESEARCH);

    double step_length;

    if(lineserach == SLEQP_LINESEARCH_EXACT)
    {
      SLEQP_CALL(sleqp_linesearch_trial_step_exact(solver->linesearch,
                                                   solver->cauchy_step,
                                                   solver->cauchy_hessian_step,
                                                   cauchy_merit_value,
                                                   solver->newton_step,
                                                   solver->newton_hessian_step,
                                                   solver->multipliers,
                                                   solver->trial_step,
                                                   &step_length,
                                                   trial_merit_value));
    }
    else
    {
      assert(lineserach == SLEQP_LINESEARCH_APPROX);

      SLEQP_CALL(sleqp_linesearch_trial_step(solver->linesearch,
                                             solver->cauchy_step,
                                             solver->cauchy_hessian_step,
                                             cauchy_merit_value,
                                             solver->newton_step,
                                             solver->newton_hessian_step,
                                             solver->multipliers,
                                             solver->trial_step,
                                             &step_length,
                                             trial_merit_value));
    }
  }

#if !defined(NDEBUG)

  {
    double actual_quadratic_merit_value;

    double func_dual = 1.;

    SLEQP_CALL(sleqp_merit_quadratic(solver->merit_data,
                                     iterate,
                                     &func_dual,
                                     solver->trial_step,
                                     solver->multipliers,
                                     solver->penalty_parameter,
                                     &actual_quadratic_merit_value));

    sleqp_assert_is_eq(*trial_merit_value,
                       actual_quadratic_merit_value,
                       eps);
  }

#endif

  SLEQP_CALL(sleqp_sparse_vector_add(sleqp_iterate_get_primal(iterate),
                                     solver->trial_step,
                                     zero_eps,
                                     solver->initial_trial_point));

  SLEQP_CALL(sleqp_sparse_vector_clip(solver->initial_trial_point,
                                      sleqp_problem_var_lb(problem),
                                      sleqp_problem_var_ub(problem),
                                      zero_eps,
                                      sleqp_iterate_get_primal(solver->trial_iterate)));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE compute_trial_point_soc(SleqpSolver* solver)
{
  const double zero_eps = sleqp_params_get(solver->params,
                                           SLEQP_PARAM_ZERO_EPS);

  SleqpProblem* problem = solver->problem;

  SleqpIterate* iterate = solver->iterate;
  SleqpIterate* trial_iterate = solver->trial_iterate;

  SleqpSparseVec* trial_point = sleqp_iterate_get_primal(trial_iterate);

  SLEQP_CALL(sleqp_soc_compute_trial_point(solver->soc_data,
                                           solver->aug_jacobian,
                                           iterate,
                                           solver->trial_step,
                                           trial_iterate,
                                           solver->soc_trial_point,
                                           &solver->soc_step_norm));

  SLEQP_CALL(sleqp_sparse_vector_clip(solver->soc_trial_point,
                                      sleqp_problem_var_lb(problem),
                                      sleqp_problem_var_ub(problem),
                                      zero_eps,
                                      trial_point));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE set_func_value(SleqpSolver* solver,
                                    SleqpIterate* iterate,
                                    SLEQP_VALUE_REASON reason)
{
  SleqpProblem* problem = solver->problem;

  int func_grad_nnz = 0;
  int cons_val_nnz = 0;
  int cons_jac_nnz = 0;

  SLEQP_CALL(sleqp_problem_set_value(problem,
                                     sleqp_iterate_get_primal(iterate),
                                     reason,
                                     &func_grad_nnz,
                                     &cons_val_nnz,
                                     &cons_jac_nnz));

  SLEQP_CALL(sleqp_sparse_vector_reserve(sleqp_iterate_get_func_grad(iterate),
                                         func_grad_nnz));

  SLEQP_CALL(sleqp_sparse_vector_reserve(sleqp_iterate_get_cons_val(iterate),
                                         cons_val_nnz));

  SLEQP_CALL(sleqp_sparse_matrix_reserve(sleqp_iterate_get_cons_jac(iterate),
                                         cons_jac_nnz));

  return SLEQP_OKAY;
}

#define HEADER_FORMAT "%10s |%14s |%14s |%14s |%14s |%14s |%14s |%14s |%14s |%14s |%14s |%14s |%14s |%14s | %18s"

#define LINE_FORMAT SLEQP_FORMAT_BOLD "%10d " SLEQP_FORMAT_RESET "|%14e |%14e |%14e |%14e |%14e |%14e |%14s |%14e |%14e |%14e |%14s |%14e |%14e | %18s"

#define INITIAL_LINE_FORMAT SLEQP_FORMAT_BOLD "%10d " SLEQP_FORMAT_RESET "|%14e |%14e |%14e |%14s |%14s |%14e |%14s |%14s |%14s |%14s |%14s |%14s |%14s | %18s"

static SLEQP_RETCODE print_header()
{
  sleqp_log_info(HEADER_FORMAT,
                 "Iteration",
                 "Func val",
                 "Merit val",
                 "Feas res",
                 "Slack res",
                 "Stat res",
                 "Penalty",
                 "Working set",
                 "LP tr",
                 "EQP tr",
                 "LP cond",
                 "Jac cond",
                 "Primal step",
                 "Dual step",
                 "Step type");

  return SLEQP_OKAY;
}

static SLEQP_RETCODE print_initial_line(SleqpSolver* solver)
{
  sleqp_log_info(INITIAL_LINE_FORMAT,
                 solver->iteration,
                 sleqp_iterate_get_func_val(solver->iterate),
                 solver->current_merit_value,
                 solver->feasibility_residuum,
                 "",
                 "",
                 solver->penalty_parameter,
                 "",
                 "",
                 "",
                 "",
                 "",
                 "",
                 "",
                 "");

  return SLEQP_OKAY;
}

static SLEQP_RETCODE print_line(SleqpSolver* solver)
{
  bool exact = false;
  double basis_condition, aug_jac_condition;

  SLEQP_CALL(sleqp_lpi_get_basis_condition(solver->lp_interface,
                                           &exact,
                                           &basis_condition));

  SLEQP_CALL(sleqp_aug_jacobian_get_condition_estimate(solver->aug_jacobian,
                                                       &aug_jac_condition));

  const char* steptype_descriptions[] = {
    [SLEQP_STEPTYPE_NONE] = "",
    [SLEQP_STEPTYPE_ACCEPTED] = "Accepted",
    [SLEQP_STEPTYPE_ACCEPTED_FULL] = "Accepted (full)",
    [SLEQP_STEPTYPE_ACCEPTED_SOC] = "Accepted SOC",
    [SLEQP_STEPTYPE_REJECTED] = "Rejected"
  };

  char jac_condition_buf[1024];

  if(aug_jac_condition != SLEQP_NONE)
  {
    snprintf(jac_condition_buf,
             1024,
             "%14e",
             aug_jac_condition);
  }
  else
  {
    snprintf(jac_condition_buf,
             1024,
             "-");
  }

  char working_set_buf[1024];

  SleqpWorkingSet* working_set = sleqp_iterate_get_working_set(solver->iterate);

  SleqpWorkingSet* trial_working_set = sleqp_iterate_get_working_set(solver->trial_iterate);

  if(sleqp_working_set_eq(working_set, trial_working_set))
  {
    snprintf(working_set_buf,
             1024,
             "--");
  }
  else
  {
    snprintf(working_set_buf,
             1024,
             "%dv/%dc",
             sleqp_working_set_num_active_vars(working_set),
             sleqp_working_set_num_active_cons(working_set));
  }

  sleqp_log_info(LINE_FORMAT,
                 solver->iteration,
                 sleqp_iterate_get_func_val(solver->iterate),
                 solver->current_merit_value,
                 solver->feasibility_residuum,
                 solver->slackness_residuum,
                 solver->stationarity_residuum,
                 solver->penalty_parameter,
                 working_set_buf,
                 solver->lp_trust_radius,
                 solver->trust_radius,
                 basis_condition,
                 jac_condition_buf,
                 solver->primal_diff_norm,
                 solver->dual_diff_norm,
                 steptype_descriptions[solver->last_step_type]);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE compute_step_lengths(SleqpSolver* solver)
{
  SleqpIterate* iterate = solver->iterate;
  SleqpIterate* trial_iterate = solver->trial_iterate;

  const double zero_eps = sleqp_params_get(solver->params, SLEQP_PARAM_ZERO_EPS);

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(sleqp_iterate_get_primal(iterate),
                                            sleqp_iterate_get_primal(trial_iterate),
                                            1.,
                                            -1,
                                            zero_eps,
                                            solver->primal_diff));

  solver->primal_diff_norm = sleqp_sparse_vector_norm(solver->primal_diff);

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(sleqp_iterate_get_cons_dual(iterate),
                                            sleqp_iterate_get_cons_dual(trial_iterate),
                                            1.,
                                            -1,
                                            zero_eps,
                                            solver->cons_dual_diff));

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(sleqp_iterate_get_vars_dual(iterate),
                                            sleqp_iterate_get_vars_dual(trial_iterate),
                                            1.,
                                            -1,
                                            zero_eps,
                                            solver->vars_dual_diff));

  solver->dual_diff_norm = 0.;

  solver->dual_diff_norm += sleqp_sparse_vector_norm_sq(solver->cons_dual_diff);
  solver->dual_diff_norm += sleqp_sparse_vector_norm_sq(solver->vars_dual_diff);

  solver->dual_diff_norm = sqrt(solver->dual_diff_norm);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE check_derivative(SleqpSolver* solver)
{
  SleqpOptions* options = solver->options;
  SleqpIterate* iterate = solver->iterate;

  const SLEQP_DERIV_CHECK deriv_check = sleqp_options_get_int(options,
                                                              SLEQP_OPTION_INT_DERIV_CHECK);

  if(deriv_check & SLEQP_DERIV_CHECK_FIRST)
  {
    SLEQP_CALL(sleqp_deriv_check_first_order(solver->deriv_check, iterate));
  }

  if(deriv_check & SLEQP_DERIV_CHECK_SECOND_EXHAUSTIVE)
  {
    SLEQP_CALL(sleqp_deriv_check_second_order_exhaustive(solver->deriv_check, iterate));
  }
  else if(deriv_check & SLEQP_DERIV_CHECK_SECOND_SIMPLE)
  {
    SLEQP_CALL(sleqp_deriv_check_second_order_simple(solver->deriv_check, iterate));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE reject_step(SleqpSolver* solver)
{
  SleqpIterate* iterate = solver->iterate;

  SLEQP_CALL(set_func_value(solver,
                            iterate,
                            SLEQP_VALUE_REASON_REJECTED_ITERATE));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE accept_step(SleqpSolver* solver)
{
  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;
  SleqpIterate* trial_iterate = solver->trial_iterate;

  SLEQP_CALL(set_func_value(solver,
                            trial_iterate,
                            SLEQP_VALUE_REASON_ACCEPTED_ITERATE));

  // get the remaining data to fill the iterate

  SLEQP_CALL(sleqp_problem_eval(problem,
                                NULL,
                                NULL,
                                sleqp_iterate_get_func_grad(trial_iterate),
                                sleqp_iterate_get_cons_val(trial_iterate),
                                sleqp_iterate_get_cons_jac(trial_iterate)));

  if(solver->bfgs_data)
  {
    SLEQP_CALL(sleqp_bfgs_push(solver->bfgs_data,
                               solver->iterate,
                               solver->trial_iterate,
                               solver->multipliers));
  }

  if(solver->sr1_data)
  {
    SLEQP_CALL(sleqp_sr1_push(solver->sr1_data,
                              solver->iterate,
                              solver->trial_iterate,
                              solver->multipliers));
  }

  // perform simple swaps
  solver->trial_iterate = iterate;
  solver->iterate = trial_iterate;

  // ensure that the unscaled iterate is kept up to date
  if(solver->scaling_data || solver->preprocessor)
  {
    solver->restore_original_iterate = true;
  }
  else
  {
    solver->original_iterate = solver->iterate;
  }

  SleqpCallbackHandler* handler = solver->callback_handlers[SLEQP_SOLVER_EVENT_ACCEPTED_ITERATE];

  if(sleqp_callback_handler_size(handler) != 0)
  {
    SLEQP_CALL(restore_original_iterate(solver));

    SLEQP_CALLBACK_HANDLER_EXECUTE(handler,
                                   SLEQP_ACCEPTED_ITERATE,
                                   solver,
                                   solver->original_iterate);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE perform_iteration(SleqpSolver* solver,
                                       bool* optimal)
{
  *optimal = false;

  SLEQP_CALL(sleqp_lpi_set_time_limit(solver->lp_interface, remaining_time(solver)));

  const SleqpOptions* options = solver->options;

  const bool quadratic_model = sleqp_options_get_bool(options,
                                                      SLEQP_OPTION_BOOL_USE_QUADRATIC_MODEL);

  const bool perform_newton_step = quadratic_model &&
    sleqp_options_get_bool(options, SLEQP_OPTION_BOOL_PERFORM_NEWTON_STEP);

  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;
  SleqpIterate* trial_iterate = solver->trial_iterate;

  const int num_constraints = sleqp_problem_num_constraints(problem);

  assert(sleqp_sparse_vector_is_boxed(sleqp_iterate_get_primal(iterate),
                                      sleqp_problem_var_lb(problem),
                                      sleqp_problem_var_ub(problem)));

  const double eps = sleqp_params_get(solver->params,
                                      SLEQP_PARAM_EPS);

  const double zero_eps = sleqp_params_get(solver->params,
                                           SLEQP_PARAM_ZERO_EPS);

  const double accepted_reduction = sleqp_params_get(solver->params,
                                                     SLEQP_PARAM_ACCEPTED_REDUCTION);

  double exact_iterate_value, model_iterate_value;

  {
    SLEQP_CALL(sleqp_merit_func(solver->merit_data,
                                iterate,
                                solver->penalty_parameter,
                                &exact_iterate_value));

    solver->current_merit_value = model_iterate_value = exact_iterate_value;
  }

  SLEQP_CALL(set_residuum(solver));

  if(solver->iteration == 0)
  {
    SLEQP_CALL(print_initial_line(solver));
  }

  double model_trial_value;

  bool full_step;

  // Derivative check
  SLEQP_CALL(check_derivative(solver));

  // Optimality check with respect to scaled problem
  {
    if(sleqp_iterate_is_optimal(iterate,
                                solver->params,
                                solver->feasibility_residuum,
                                solver->slackness_residuum,
                                solver->stationarity_residuum))
    {
      *optimal = true;
    }
  }

  if(*optimal)
  {
    return SLEQP_OKAY;
  }

  // Step computation
  if(perform_newton_step)
  {
    SLEQP_CALL(compute_trial_point_newton(solver,
                                          &model_trial_value,
                                          &full_step));
  }
  else
  {
    SLEQP_CALL(compute_trial_point_simple(solver,
                                          &model_trial_value,
                                          quadratic_model,
                                          &full_step));
  }

  SLEQP_CALL(compute_step_lengths(solver));

  double model_reduction = model_iterate_value - model_trial_value;

  sleqp_assert_is_geq(model_reduction, 0., zero_eps);

  SLEQP_CALL(set_func_value(solver, trial_iterate, SLEQP_VALUE_REASON_TRYING_ITERATE));

  double func_val;

  SLEQP_CALL(sleqp_problem_eval(problem,
                                NULL,
                                &func_val,
                                NULL,
                                sleqp_iterate_get_cons_val(trial_iterate),
                                NULL));

  SLEQP_CALL(sleqp_iterate_set_func_val(trial_iterate, func_val));

  double actual_reduction = 0.;

  {
    double exact_trial_value;

    SLEQP_CALL(sleqp_merit_func(solver->merit_data,
                                trial_iterate,
                                solver->penalty_parameter,
                                &exact_trial_value));

    actual_reduction = exact_iterate_value - exact_trial_value;

    sleqp_log_debug("Current merit function value: %e, trial merit function value: %e",
                    exact_iterate_value,
                    exact_trial_value);

  }

  double reduction_ratio = 1.;

  if(actual_reduction != model_reduction)
  {
    reduction_ratio = actual_reduction / model_reduction;
  }

  sleqp_log_debug("Reduction ratio: %e, actual: %e, predicted: %e",
                  reduction_ratio,
                  actual_reduction,
                  model_reduction);

  const double trial_step_infnorm = sleqp_sparse_vector_inf_norm(solver->trial_step);
  const double cauchy_step_infnorm = sleqp_sparse_vector_inf_norm(solver->cauchy_step);

  double trial_step_norm = sleqp_sparse_vector_norm(solver->trial_step);

  sleqp_log_debug("Trial step norm: %e", trial_step_norm);

  solver->boundary_step = sleqp_is_geq(trial_step_norm,
                                       solver->trust_radius,
                                       eps);

  bool step_accepted = true;

  solver->last_step_type = SLEQP_STEPTYPE_REJECTED;

  if(reduction_ratio >= accepted_reduction)
  {
    sleqp_log_debug("Trial step accepted");

    if(full_step)
    {
      solver->last_step_type = SLEQP_STEPTYPE_ACCEPTED_FULL;
    }
    else
    {
      solver->last_step_type = SLEQP_STEPTYPE_ACCEPTED;
    }
  }
  else
  {
    sleqp_log_debug("Trial step rejected");

    step_accepted = false;

    const bool perform_soc = sleqp_options_get_bool(options,
                                                    SLEQP_OPTION_BOOL_PERFORM_SOC);

    if((num_constraints > 0) && perform_soc)
    {
      sleqp_log_debug("Computing second-order correction");

      SLEQP_CALL(compute_trial_point_soc(solver));

      solver->boundary_step = sleqp_is_geq(solver->soc_step_norm,
                                           solver->trust_radius,
                                           eps);

      SLEQP_CALL(set_func_value(solver,
                                trial_iterate,
                                SLEQP_VALUE_REASON_TRYING_SOC_ITERATE));

      double func_val;

      SLEQP_CALL(sleqp_problem_eval(problem,
                                    NULL,
                                    &func_val,
                                    NULL,
                                    sleqp_iterate_get_cons_val(trial_iterate),
                                    NULL));

      SLEQP_CALL(sleqp_iterate_set_func_val(trial_iterate, func_val));

      double soc_exact_trial_value;

      SLEQP_CALL(sleqp_merit_func(solver->merit_data,
                                  trial_iterate,
                                  solver->penalty_parameter,
                                  &soc_exact_trial_value));

      const double soc_actual_reduction = exact_iterate_value - soc_exact_trial_value;

      reduction_ratio = 1.;

      // The denominator of the SOC reduction ratio is the quadratic reduction
      // with respect to the original (not the corrected) trial step
      if(soc_actual_reduction != model_reduction)
      {
        reduction_ratio = soc_actual_reduction / model_reduction;
      }

      sleqp_log_debug("SOC Reduction ratio: %e, actual: %e, predicted: %e",
                      reduction_ratio,
                      soc_actual_reduction,
                      model_reduction);

      if(reduction_ratio >= accepted_reduction)
      {
        step_accepted = true;
      }

      if(step_accepted)
      {
        solver->last_step_type = SLEQP_STEPTYPE_ACCEPTED_SOC;
        sleqp_log_debug("Second-order correction accepted");
      }
      else
      {
        sleqp_log_debug("Second-order correction rejected");
      }
    }
  }

  ++solver->iteration;

  if(solver->iteration % 25 == 0)
  {
    SLEQP_CALL(print_header());
  }

  SLEQP_CALL(print_line(solver));

  // update trust radii, penalty parameter
  {
    if(perform_newton_step)
    {
      SLEQP_CALL(update_trust_radius(solver,
                                     reduction_ratio,
                                     step_accepted,
                                     trial_step_norm));
    }

    SLEQP_CALL(update_lp_trust_radius(step_accepted,
                                      trial_step_infnorm,
                                      cauchy_step_infnorm,
                                      solver->cauchy_step_length,
                                      zero_eps,
                                      &(solver->lp_trust_radius)));

    SLEQP_CALL(sleqp_update_penalty(problem,
                                    iterate,
                                    solver->cauchy_data,
                                    &(solver->penalty_parameter)));
  }

  // update current iterate

  if(step_accepted)
  {
    SLEQP_CALL(accept_step(solver));
  }
  else
  {
    SLEQP_CALL(reject_step(solver));
  }

  SLEQP_CALLBACK_HANDLER_EXECUTE(solver->callback_handlers[SLEQP_SOLVER_EVENT_PERFORMED_ITERATION],
                                 SLEQP_PERFORMED_ITERATION,
                                 solver);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE solver_print_timer(SleqpTimer* timer,
                                        const char* message,
                                        double total_elapsed)
{
  const int buf_size = 4096;
  char buffer[buf_size];

  const int num_runs = sleqp_timer_get_num_runs(timer);
  const double avg_time = sleqp_timer_get_avg(timer);
  const double total_time = sleqp_timer_get_ttl(timer);
  const double percent = (total_time / total_elapsed) * 100.;

  snprintf(buffer,
           buf_size,
           "%30s: %5d (%.6fs avg, %8.2fs total = %5.2f%%)",
           message,
           num_runs,
           avg_time,
           total_time,
           percent);

  sleqp_log_info(buffer);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE solver_print_stats(SleqpSolver* solver,
                                        double violation)
{
  const char* descriptions[] = {
    [SLEQP_FEASIBLE]   = SLEQP_FORMAT_BOLD SLEQP_FORMAT_YELLOW "feasible"   SLEQP_FORMAT_RESET,
    [SLEQP_OPTIMAL]    = SLEQP_FORMAT_BOLD SLEQP_FORMAT_GREEN  "optimal"    SLEQP_FORMAT_RESET,
    [SLEQP_INFEASIBLE] = SLEQP_FORMAT_BOLD SLEQP_FORMAT_RED    "infeasible" SLEQP_FORMAT_RESET,
    [SLEQP_INVALID]    = SLEQP_FORMAT_BOLD SLEQP_FORMAT_RED    "Invalid"    SLEQP_FORMAT_RESET
  };

  SleqpFunc* original_func = sleqp_problem_func(solver->original_problem);
  SleqpFunc* func = sleqp_problem_func(solver->problem);

  const bool with_hessian = !(solver->sr1_data || solver->bfgs_data);

  sleqp_log_info(SLEQP_FORMAT_BOLD "%30s: %s" SLEQP_FORMAT_RESET,
                 "Solution status",
                 descriptions[solver->status]);

  if(solver->scaling_data)
  {
    double unscaled_violation;

    SLEQP_CALL(sleqp_iterate_feasibility_residuum(solver->original_problem,
                                                  solver->original_iterate,
                                                  &unscaled_violation));

    sleqp_log_info(SLEQP_FORMAT_BOLD "%30s:     %5.10e" SLEQP_FORMAT_RESET,
                   "Scaled objective value",
                   sleqp_iterate_get_func_val(solver->iterate));

    sleqp_log_info(SLEQP_FORMAT_BOLD "%30s:     %5.10e" SLEQP_FORMAT_RESET,
                   "Scaled violation",
                   violation);

    sleqp_log_info("%30s:     %5.10e",
                   "Original objective value",
                   sleqp_iterate_get_func_val(solver->original_iterate));

    sleqp_log_info("%30s:     %5.10e",
                   "Original violation",
                   unscaled_violation);

  }
  else
  {
    sleqp_log_info(SLEQP_FORMAT_BOLD "%30s:     %5.10e" SLEQP_FORMAT_RESET,
                   "Objective value",
                   sleqp_iterate_get_func_val(solver->iterate));

    sleqp_log_info(SLEQP_FORMAT_BOLD "%30s:     %5.10e" SLEQP_FORMAT_RESET,
                   "Violation",
                   violation);
  }

  sleqp_log_info("%30s: %5d",
                 "Iterations",
                 solver->iteration);

  SLEQP_CALL(solver_print_timer(sleqp_func_get_set_timer(original_func),
                                "Setting function values",
                                solver->elapsed_seconds));

  SLEQP_CALL(solver_print_timer(sleqp_func_get_val_timer(original_func),
                                "Function evaluations",
                                solver->elapsed_seconds));

  SLEQP_CALL(solver_print_timer(sleqp_func_get_grad_timer(original_func),
                                "Gradient evaluations",
                                solver->elapsed_seconds));

  SLEQP_CALL(solver_print_timer(sleqp_func_get_cons_val_timer(original_func),
                                "Constraint evaluations",
                                solver->elapsed_seconds));

  SLEQP_CALL(solver_print_timer(sleqp_func_get_cons_jac_timer(original_func),
                                "Jacobian evaluations",
                                solver->elapsed_seconds));

  if(with_hessian)
  {
    SLEQP_CALL(solver_print_timer(sleqp_func_get_hess_timer(original_func),
                                  "Hessian products",
                                  solver->elapsed_seconds));
  }

  if(solver->bfgs_data)
  {
    SLEQP_CALL(solver_print_timer(sleqp_func_get_hess_timer(func),
                                  "BFGS products",
                                  solver->elapsed_seconds));

    SLEQP_CALL(solver_print_timer(sleqp_bfgs_update_timer(solver->bfgs_data),
                                  "BFGS updates",
                                  solver->elapsed_seconds));
  }

  if(solver->sr1_data)
  {
    SLEQP_CALL(solver_print_timer(sleqp_func_get_hess_timer(func),
                                  "SR1 products",
                                  solver->elapsed_seconds));

    SLEQP_CALL(solver_print_timer(sleqp_sr1_update_timer(solver->sr1_data),
                                  "SR1 updates",
                                  solver->elapsed_seconds));
  }

  SLEQP_CALL(solver_print_timer(sleqp_aug_jacobian_get_factorization_timer(solver->aug_jacobian),
                                "Factorizations",
                                solver->elapsed_seconds));

  SLEQP_CALL(solver_print_timer(sleqp_aug_jacobian_get_substitution_timer(solver->aug_jacobian),
                                "Substitutions",
                                solver->elapsed_seconds));

  SLEQP_CALL(solver_print_timer(sleqp_lpi_get_solve_timer(solver->lp_interface),
                                "Solved LPs",
                                solver->elapsed_seconds));

  SLEQP_CALL(solver_print_timer(sleqp_newton_get_timer(solver->newton_data),
                                "Solved EQPs",
                                solver->elapsed_seconds));

  SLEQP_CALL(solver_print_timer(sleqp_linesearch_get_timer(solver->linesearch),
                                "Line searches",
                                solver->elapsed_seconds));

  sleqp_log_info("%30s: %8.2fs", "Solving time", solver->elapsed_seconds);

  if(solver->status == SLEQP_INFEASIBLE)
  {
    SLEQP_CALL(restore_original_iterate(solver));

    SLEQP_CALL(sleqp_violation_values(solver->original_problem,
                                      solver->original_iterate,
                                      solver->original_violation));

    sleqp_log_info("Violations with respect to original problem: ");

    for(int index = 0; index < solver->original_violation->nnz; ++index)
    {
      sleqp_log_info("(%d) = %e",
                     solver->original_violation->indices[index],
                     solver->original_violation->data[index]);
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solver_solve(SleqpSolver* solver,
                                 int max_num_iterations,
                                 double time_limit)
{
  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;

  const int num_variables = sleqp_problem_num_variables(problem);
  const int num_constraints = sleqp_problem_num_constraints(problem);

  const double feas_eps = sleqp_params_get(solver->params,
                                           SLEQP_PARAM_FEASIBILITY_TOL);

  sleqp_log_info("Solving a problem with %d variables, %d constraints",
                 num_variables,
                 num_constraints);

  SLEQP_CALL(sleqp_set_and_evaluate(problem,
                                    iterate,
                                    SLEQP_VALUE_REASON_INIT));

  // ensure that the unscaled iterate is initialized
  if(solver->scaling_data)
  {
    SLEQP_CALL(sleqp_iterate_copy(iterate,
                                  solver->original_iterate));

    SLEQP_CALL(sleqp_unscale_iterate(solver->scaling_data,
                                     solver->original_iterate));
  }

  // Warnings
  {
    const SLEQP_DERIV_CHECK deriv_check = sleqp_options_get_int(solver->options,
                                                                SLEQP_OPTION_INT_DERIV_CHECK);

    const bool inexact_hessian = (solver->bfgs_data || solver->sr1_data);

    const int hessian_check_flags = SLEQP_DERIV_CHECK_SECOND_EXHAUSTIVE | SLEQP_DERIV_CHECK_SECOND_SIMPLE;

    const bool hessian_check = (deriv_check & hessian_check_flags);

    if(inexact_hessian && hessian_check)
    {
      sleqp_log_warn("Enabled second order derivative check while using a quasi-Newton method");
    }
  }

  {
    double total_violation;

    SLEQP_CALL(sleqp_violation_one_norm(problem,
                                        sleqp_iterate_get_cons_val(iterate),
                                        &total_violation));

    const double func_val = sleqp_iterate_get_func_val(iterate);

    if(total_violation >= 10. * func_val)
    {
      sleqp_log_warn("Problem is badly scaled, constraint violation %g significantly exceeds function value of %g",
                     total_violation,
                     func_val);
    }
  }

  solver->status = SLEQP_INVALID;

  SLEQP_CALL(sleqp_timer_reset(solver->elapsed_timer));

  solver->time_limit = time_limit;
  solver->abort_next = false;

  solver->iteration = 0;
  solver->elapsed_seconds = 0.;
  solver->last_step_type = SLEQP_STEPTYPE_NONE;

  const double deadpoint_bound = sleqp_params_get(solver->params,
                                                  SLEQP_PARAM_DEADPOINT_BOUND);

  bool reached_deadpoint = false;

  SLEQP_CALL(print_header());

  // main solving loop
  while(true)
  {
    if(time_limit != SLEQP_NONE)
    {
      if(solver->elapsed_seconds >= time_limit)
      {
        sleqp_log_info("Exhausted time limit, terminating");
        break;
      }
    }

    if(max_num_iterations != SLEQP_NONE &&
       solver->iteration >= max_num_iterations)
    {
      sleqp_log_info("Reached iteration limit, terminating");
      break;
    }

    if(solver->abort_next)
    {
      sleqp_log_info("Abortion requested, terminating");
      break;
    }

    bool optimal;

    SLEQP_CALL(sleqp_timer_start(solver->elapsed_timer));

    SLEQP_CALL(perform_iteration(solver, &optimal));

    SLEQP_CALL(sleqp_timer_stop(solver->elapsed_timer));

    solver->elapsed_seconds = sleqp_timer_get_ttl(solver->elapsed_timer);

    if(solver->lp_trust_radius <= deadpoint_bound ||
       solver->trust_radius <= deadpoint_bound)
    {
      reached_deadpoint = true;
      break;
    }

    if(optimal)
    {
      sleqp_log_debug("Achieved optimality");
      solver->status = SLEQP_OPTIMAL;
      break;
    }
  }

  if(reached_deadpoint)
  {
    sleqp_log_warn("Reached dead point");
  }

  double violation;

  SLEQP_CALL(sleqp_iterate_feasibility_residuum(solver->problem,
                                                solver->iterate,
                                                &violation));

  if(solver->status != SLEQP_OPTIMAL)
  {
    const bool feasible = sleqp_iterate_is_feasible(iterate,
                                                    solver->feasibility_residuum,
                                                    feas_eps);

    solver->status = feasible ? SLEQP_FEASIBLE : SLEQP_INFEASIBLE;
  }

  SLEQP_CALLBACK_HANDLER_EXECUTE(solver->callback_handlers[SLEQP_SOLVER_EVENT_FINISHED],
                                 SLEQP_FINISHED,
                                 solver,
                                 solver->original_iterate);

  SLEQP_CALL(solver_print_stats(solver, violation));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solver_get_real_state(const SleqpSolver* solver,
                                          SLEQP_SOLVER_STATE_REAL state,
                                          double* value)
{
  double min_rayleigh, max_rayleigh;

  SLEQP_CALL(sleqp_newton_current_rayleigh(solver->newton_data,
                                           &min_rayleigh,
                                           &max_rayleigh));

  switch(state)
  {
  case SLEQP_SOLVER_STATE_REAL_TRUST_RADIUS:
    (*value) = solver->trust_radius;
    break;
  case SLEQP_SOLVER_STATE_REAL_LP_TRUST_RADIUS:
    (*value) = solver->lp_trust_radius;
    break;
  case SLEQP_SOLVER_STATE_REAL_SCALED_FUNC_VAL:
    (*value) = sleqp_iterate_get_func_val(solver->iterate);
    break;
  case SLEQP_SOLVER_STATE_REAL_SCALED_MERIT_VAL:
    (*value) = solver->current_merit_value;
    break;
  case SLEQP_SOLVER_STATE_REAL_SCALED_FEAS_RES:
    (*value) = solver->feasibility_residuum;
    break;
  case SLEQP_SOLVER_STATE_REAL_SCALED_STAT_RES:
    (*value) = solver->stationarity_residuum;
    break;
  case SLEQP_SOLVER_STATE_REAL_SCALED_SLACK_RES:
    (*value) = solver->slackness_residuum;
    break;
  case SLEQP_SOLVER_STATE_REAL_PENALTY_PARAM:
    (*value) = solver->penalty_parameter;
    break;
  case SLEQP_SOLVER_STATE_REAL_MIN_RAYLEIGH:
    (*value) = min_rayleigh;
    break;
  case SLEQP_SOLVER_STATE_REAL_MAX_RAYLEIGH:
    (*value) = max_rayleigh;
    break;
  default:
    sleqp_log_error("Invalid state requested (%d)", value);
    return SLEQP_ILLEGAL_ARGUMENT;
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solver_get_int_state(const SleqpSolver* solver,
                                         SLEQP_SOLVER_STATE_INT state,
                                         int* value)
{
  switch(state)
  {
  case SLEQP_SOLVER_STATE_INT_LAST_STEP_TYPE:
    (*value) = solver->last_step_type;
    break;
  case SLEQP_SOLVER_STATE_INT_LAST_STEP_ON_BDRY:
    (*value) = solver->boundary_step;
    break;
  default:
    sleqp_log_error("Invalid state requested (%d)", value);
    return SLEQP_ILLEGAL_ARGUMENT;
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solver_get_vec_state(const SleqpSolver* solver,
                                         SLEQP_SOLVER_STATE_VEC value,
                                         SleqpSparseVec* result)
{
  const double zero_eps = sleqp_params_get(solver->params,
                                           SLEQP_PARAM_ZERO_EPS);

  switch(value)
  {
  case SLEQP_SOLVER_STATE_VEC_SCALED_STAT_RESIDUALS:
    SLEQP_CALL(sleqp_iterate_stationarity_residuals(solver->problem,
                                                    solver->iterate,
                                                    solver->dense_cache,
                                                    result,
                                                    zero_eps));
    break;
  case SLEQP_SOLVER_STATE_VEC_SCALED_FEAS_RESIDUALS:
    SLEQP_CALL(sleqp_violation_values(solver->problem,
                                      solver->iterate,
                                      result));
    break;
  case SLEQP_SOLVER_STATE_VEC_SCALED_CONS_SLACK_RESIDUALS:
    SLEQP_CALL(sleqp_iterate_cons_slackness_residuals(solver->problem,
                                                      solver->iterate,
                                                      result,
                                                      zero_eps));
    break;
  case SLEQP_SOLVER_STATE_VEC_SCALED_VAR_SLACK_RESIDUALS:
    SLEQP_CALL(sleqp_iterate_vars_slackness_residuals(solver->problem,
                                                      solver->iterate,
                                                      result,
                                                      zero_eps));
    break;
  default:
    sleqp_log_error("Invalid state requested (%d)", value);
    return SLEQP_ILLEGAL_ARGUMENT;
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solver_get_solution(SleqpSolver* solver,
                                        SleqpIterate** iterate)
{
  SLEQP_CALL(restore_original_iterate(solver));

  (*iterate) = solver->original_iterate;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solver_get_violated_constraints(SleqpSolver* solver,
                                                    SleqpIterate* iterate,
                                                    int* violated_constraints,
                                                    int* num_violated_constraints)
{
  const double feas_eps = sleqp_params_get(solver->params,
                                           SLEQP_PARAM_FEASIBILITY_TOL);

  SLEQP_CALL(sleqp_iterate_get_violated_constraints(solver->original_problem,
                                                    iterate,
                                                    violated_constraints,
                                                    num_violated_constraints,
                                                    feas_eps));

  return SLEQP_OKAY;
}

SLEQP_STATUS sleqp_solver_get_status(const SleqpSolver* solver)
{
  return solver->status;
}

SLEQP_RETCODE sleqp_solver_reset(SleqpSolver* solver)
{
  const int num_variables = sleqp_problem_num_variables(solver->problem);

  // initial trust region radii as suggested,
  // penalty parameter as suggested:

  solver->trust_radius = 1.;
  solver->lp_trust_radius = .8 * (solver->trust_radius) * sqrt((double) num_variables);

  solver->penalty_parameter = 10.;

  if(solver->bfgs_data)
  {
    SLEQP_CALL(sleqp_bfgs_reset(solver->bfgs_data));
  }

  if(solver->sr1_data)
  {
    SLEQP_CALL(sleqp_sr1_reset(solver->sr1_data));
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solver_abort(SleqpSolver* solver)
{
  solver->abort_next = true;

  return SLEQP_OKAY;
}

int sleqp_solver_get_iterations(const SleqpSolver* solver)
{
  return solver->iteration;
}

double sleqp_solver_get_elapsed_seconds(const SleqpSolver* solver)
{
  return solver->elapsed_seconds;
}

SLEQP_RETCODE sleqp_solver_add_callback(SleqpSolver* solver,
                                        SLEQP_SOLVER_EVENT solver_event,
                                        void* callback_func,
                                        void* callback_data)
{
  if(solver_event < 0 || solver_event >= SLEQP_SOLVER_NUM_EVENTS)
  {
    sleqp_log_error("Invalid callback");
    return SLEQP_ILLEGAL_ARGUMENT;
  }

  SLEQP_CALL(sleqp_callback_handler_add(solver->callback_handlers[solver_event],
                                        callback_func,
                                        callback_data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solver_remove_callback(SleqpSolver* solver,
                                           SLEQP_SOLVER_EVENT solver_event,
                                           void* callback_func,
                                           void* callback_data)
{
  if(solver_event < 0 || solver_event >= SLEQP_SOLVER_NUM_EVENTS)
  {
    sleqp_log_error("Invalid callback");
    return SLEQP_ILLEGAL_ARGUMENT;
  }

  SLEQP_CALL(sleqp_callback_handler_remove(solver->callback_handlers[solver_event],
                                           callback_func,
                                           callback_data));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE solver_free(SleqpSolver** star)
{
  SleqpSolver* solver = *star;

  if(!solver)
  {
    return SLEQP_OKAY;
  }

  sleqp_free(&solver->dense_cache);

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->soc_trial_point));

  SLEQP_CALL(sleqp_soc_data_release(&solver->soc_data));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->vars_dual_diff));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->cons_dual_diff));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->primal_diff));

  if(solver->callback_handlers)
  {
    for(int i = 0; i < SLEQP_SOLVER_NUM_EVENTS; ++i)
    {
      SLEQP_CALL(sleqp_callback_handler_release(solver->callback_handlers + i));
    }
  }

  sleqp_free(&solver->callback_handlers);

  SLEQP_CALL(sleqp_working_set_release(&solver->parametric_original_working_set));

  SLEQP_CALL(sleqp_parametric_solver_release(&solver->parametric_solver));

  SLEQP_CALL(sleqp_linesearch_release(&solver->linesearch));

  SLEQP_CALL(sleqp_merit_data_release(&solver->merit_data));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->estimation_residuals));

  SLEQP_CALL(sleqp_dual_estimation_free(&solver->estimation_data));

  SLEQP_CALL(sleqp_aug_jacobian_release(&solver->aug_jacobian));

  SLEQP_CALL(sleqp_sparse_factorization_release(&solver->factorization));

  SLEQP_CALL(sleqp_iterate_release(&solver->trial_iterate));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->initial_trial_point));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->trial_step));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->newton_hessian_step));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->newton_step));

  SLEQP_CALL(sleqp_newton_data_release(&solver->newton_data));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->multipliers));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->cauchy_hessian_step));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->cauchy_step));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->cauchy_direction));

  SLEQP_CALL(sleqp_cauchy_release(&solver->cauchy_data));

  SLEQP_CALL(sleqp_lpi_free(&solver->lp_interface));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->original_violation));

  SLEQP_CALL(sleqp_iterate_release(&solver->iterate));

  SLEQP_CALL(sleqp_deriv_checker_free(&solver->deriv_check));

  SLEQP_CALL(sleqp_options_release(&solver->options));
  SLEQP_CALL(sleqp_params_release(&solver->params));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->primal));
  SLEQP_CALL(sleqp_sparse_vector_free(&solver->scaled_primal));

  SLEQP_CALL(sleqp_timer_free(&solver->elapsed_timer));

  SLEQP_CALL(sleqp_iterate_release(&solver->scaled_iterate));

  if(solver->scaling_data || solver->preprocessor)
  {
    SLEQP_CALL(sleqp_iterate_release(&solver->original_iterate));
  }

  SLEQP_CALL(sleqp_problem_release(&solver->problem));

  SLEQP_CALL(sleqp_problem_scaling_release(&solver->problem_scaling));

  SLEQP_CALL(sleqp_scaling_release(&solver->scaling_data));

  SLEQP_CALL(sleqp_preprocessor_release(&solver->preprocessor));

  SLEQP_CALL(sleqp_problem_release(&solver->scaled_problem));

  SLEQP_CALL(sleqp_sr1_release(&solver->sr1_data));

  SLEQP_CALL(sleqp_bfgs_release(&solver->bfgs_data));

  SLEQP_CALL(sleqp_problem_release(&solver->original_problem));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solver_capture(SleqpSolver* solver)
{
  ++solver->refcount;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solver_release(SleqpSolver** star)
{
  SleqpSolver* solver = *star;

  if(!solver)
  {
    return SLEQP_OKAY;
  }

  if(--solver->refcount == 0)
  {
    SLEQP_CALL(solver_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
