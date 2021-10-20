#include "solver.h"

#include <assert.h>
#include <math.h>
#include <string.h>
#include <threads.h>

#include "cmp.h"
#include "defs.h"
#include "fail.h"
#include "feas.h"
#include "iterate.h"
#include "mem.h"
#include "penalty.h"
#include "scale.h"

#include "timer.h"
#include "working_step.h"
#include "util.h"

#include "lp/lpi.h"

#include "aug_jac/box_constrained_aug_jac.h"
#include "aug_jac/standard_aug_jac.h"
#include "aug_jac/unconstrained_aug_jac.h"

#include "cauchy/box_constrained_cauchy.h"
#include "cauchy/standard_cauchy.h"
#include "cauchy/unconstrained_cauchy.h"

#include "step/step_rule.h"

#define INFO_BUF_SIZE 100
#define SOLVER_INFO_BUF_SIZE 400

thread_local char lps_info[INFO_BUF_SIZE];
thread_local char fact_info[INFO_BUF_SIZE];
thread_local char solver_info[SOLVER_INFO_BUF_SIZE];

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
SLEQP_RETCODE do_restore_iterate(SleqpSolver* solver,
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

  SleqpProblem* problem = solver->problem;
  SleqpFunc* func = sleqp_problem_func(problem);
  const bool lsq = sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_LSQ;

  if(solver->scaling_data)
  {
    SLEQP_CALL(sleqp_unscale_iterate(solver->scaling_data,
                                     target,
                                     lsq));
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solver_restore_original_iterate(SleqpSolver* solver)
{
  if(solver->restore_original_iterate)
  {
    SLEQP_CALL(do_restore_iterate(solver,
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

  SLEQP_CALL(sleqp_quasi_newton_create_default(&solver->quasi_newton,
                                               func,
                                               params,
                                               options));

  if(solver->quasi_newton)
  {
    func = sleqp_quasi_newton_get_func(solver->quasi_newton);
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
        solver->status = SLEQP_STATUS_INFEASIBLE;
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
  SleqpProblem* original_problem = solver->original_problem;

  const double zero_eps = sleqp_params_get(solver->params,
                                           SLEQP_PARAM_ZERO_EPS);

  {
    SleqpSparseVec* var_lb = sleqp_problem_var_lb(original_problem);
    SleqpSparseVec* var_ub = sleqp_problem_var_ub(original_problem);

    if(!sleqp_sparse_vector_is_boxed(primal, var_lb, var_ub))
    {
      sleqp_log_warn("Initial solution violates variable bounds");
    }
  }

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
                                    original_problem,
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

static
SLEQP_RETCODE create_cauchy_solver(SleqpSolver* solver)
{
  SleqpProblem* problem = solver->problem;
  SleqpParams* params = solver->params;
  SleqpOptions* options = solver->options;

  const int num_constraints = sleqp_problem_num_constraints(problem);

  if(sleqp_problem_is_unconstrained(problem))
  {
    SLEQP_CALL(sleqp_unconstrained_cauchy_create(&solver->cauchy_data,
                                                 problem,
                                                 params));
  }
  else if(num_constraints == 0)
  {
    SLEQP_CALL(sleqp_box_constrained_cauchy_create(&solver->cauchy_data,
                                                   problem,
                                                   params));
  }
  else
  {
    const int num_variables = sleqp_problem_num_variables(problem);
    const int num_lp_variables = num_variables + 2*num_constraints;
    const int num_lp_constraints = num_constraints;

    SLEQP_CALL(sleqp_lpi_create_default_interface(&solver->lp_interface,
                                                  num_lp_variables,
                                                  num_lp_constraints,
                                                  params,
                                                  options));

    SLEQP_CALL(sleqp_standard_cauchy_create(&solver->cauchy_data,
                                            problem,
                                            params,
                                            options,
                                            solver->lp_interface));
  }

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE create_aug_jac(SleqpSolver* solver)
{
  SleqpProblem* problem = solver->problem;
  SleqpParams* params = solver->params;

  const int num_constraints = sleqp_problem_num_constraints(problem);

  if(sleqp_problem_is_unconstrained(problem))
  {
    SLEQP_CALL(sleqp_unconstrained_aug_jac_create(&solver->aug_jac,
                                                  problem));
  }
  else if(num_constraints == 0)
  {
    SLEQP_CALL(sleqp_box_constrained_aug_jac_create(&solver->aug_jac,
                                                    problem));
  }
  else
  {
    // create sparse factorization

    SLEQP_CALL(sleqp_sparse_factorization_create_default(&solver->factorization,
                                                         params));

    SLEQP_CALL(sleqp_standard_aug_jac_create(&solver->aug_jac,
                                             problem,
                                             params,
                                             solver->factorization));
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

  SLEQP_CALL(sleqp_step_rule_create_default(&solver->step_rule,
                                            solver->problem,
                                            solver->params,
                                            solver->options));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->original_violation,
                                              num_constraints));

  SLEQP_CALL(create_cauchy_solver(solver));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->cauchy_direction,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->cauchy_step,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->cauchy_hessian_step,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->multipliers,
                                              num_constraints));

  SLEQP_CALL(sleqp_working_step_create(&solver->working_step,
                                       solver->problem,
                                       params));


  SLEQP_TR_SOLVER tr_solver = sleqp_options_get_int(options,
                                                    SLEQP_OPTION_INT_TR_SOLVER);

  if(tr_solver == SLEQP_TR_SOLVER_LSQR)
  {
    SleqpFunc* func = sleqp_problem_func(problem);

    if(sleqp_func_get_type(func) != SLEQP_FUNC_TYPE_LSQ)
    {
      sleqp_log_error("LSQR solver is only available for LSQ problems");
      return SLEQP_ILLEGAL_ARGUMENT;
    }

    SLEQP_CALL(sleqp_gauss_newton_solver_create(&solver->gauss_newton_solver,
                                                solver->problem,
                                                solver->working_step,
                                                solver->params));
  }

  SLEQP_CALL(sleqp_newton_data_create(&solver->newton_data,
                                      solver->problem,
                                      solver->working_step,
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

  SLEQP_CALL(create_aug_jac(solver));

  SLEQP_CALL(sleqp_dual_estimation_create(&solver->estimation_data,
                                          solver->problem));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->estimation_residuals,
                                              num_variables));

  SLEQP_CALL(sleqp_merit_create(&solver->merit,
                                solver->problem,
                                params));

  SLEQP_CALL(sleqp_linesearch_create(&solver->linesearch,
                                     solver->problem,
                                     params,
                                     solver->merit));

  SLEQP_CALL(sleqp_polishing_create(&solver->polishing,
                                    solver->problem,
                                    solver->params));

  SLEQP_PARAMETRIC_CAUCHY parametric_cauchy = sleqp_options_get_int(solver->options,
                                                                    SLEQP_OPTION_INT_PARAMETRIC_CAUCHY);

  if(parametric_cauchy != SLEQP_PARAMETRIC_CAUCHY_DISABLED)
  {
    SLEQP_CALL(sleqp_parametric_solver_create(&solver->parametric_solver,
                                              solver->problem,
                                              solver->params,
                                              solver->options,
                                              solver->merit,
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

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&solver->soc_step,
                                              num_variables));

  SLEQP_CALL(sleqp_alloc_array(&solver->dense_cache, SLEQP_MAX(num_variables, num_constraints)));

  SLEQP_CALL(sleqp_solver_reset(solver));

  solver->time_limit = SLEQP_NONE;

  solver->abort_next = false;

  sleqp_log_debug("%s", sleqp_solver_info(solver));

  return SLEQP_OKAY;
}

static void print_solver(char* buffer,
                         int len,
                         const char* name,
                         const char* version)
{
  if(strlen(version) == 0)
  {
    snprintf(buffer, len, "%s", name);
  }
  else
  {
    snprintf(buffer, len, "%s %s", name, version);
  }
}

const char* sleqp_solver_info(const SleqpSolver* solver)
{
  if(solver->lp_interface)
  {
    print_solver(lps_info,
                 INFO_BUF_SIZE,
                 sleqp_lpi_get_name(solver->lp_interface),
                 sleqp_lpi_get_version(solver->lp_interface));
  }
  else
  {
    snprintf(lps_info,
             INFO_BUF_SIZE,
             "None");
  }

  if(solver->factorization)
  {
    print_solver(fact_info,
                 INFO_BUF_SIZE,
                 sleqp_sparse_factorization_get_name(solver->factorization),
                 sleqp_sparse_factorization_get_version(solver->factorization));
  }
  else
  {
    snprintf(fact_info,
             INFO_BUF_SIZE,
             "None");
  }

  snprintf(solver_info,
           SOLVER_INFO_BUF_SIZE,
           "Sleqp version %s [LP solver: %s] [Factorization: %s] [GitHash %s]",
           SLEQP_VERSION,
           lps_info,
           fact_info,
           SLEQP_GIT_COMMIT_HASH);

  return solver_info;
}

static SLEQP_RETCODE print_warning(SleqpSolver* solver)
{
  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;

  SleqpFunc* func = sleqp_problem_func(problem);

  SLEQP_FUNC_TYPE func_type = sleqp_func_get_type(func);

  const SLEQP_DERIV_CHECK deriv_check = sleqp_options_get_int(solver->options,
                                                              SLEQP_OPTION_INT_DERIV_CHECK);

  const int hessian_check_flags = (SLEQP_DERIV_CHECK_SECOND_EXHAUSTIVE |
                                   SLEQP_DERIV_CHECK_SECOND_SIMPLE);

  const bool hessian_check = (deriv_check & hessian_check_flags);

  if(hessian_check)
  {
    const bool inexact_hessian = (solver->quasi_newton ||
                                  func_type == SLEQP_FUNC_TYPE_LSQ);

    if (inexact_hessian)
    {
      sleqp_log_warn("Enabled second order derivative check while using a "
                     "quasi-Newton method");
    }
  }

  if(func_type == SLEQP_FUNC_TYPE_DYNAMIC)
  {
    const int check_first_flags = SLEQP_DERIV_CHECK_FIRST;

    const bool check_first = (deriv_check & check_first_flags);

    if(check_first)
    {
      sleqp_log_warn("Enabled first order derivative check while using a dynamic function");
    }
  }

  {
    double total_violation;

    SLEQP_CALL(sleqp_violation_one_norm(problem,
                                        sleqp_iterate_get_cons_val(iterate),
                                        &total_violation));

    const double func_val = sleqp_iterate_get_func_val(iterate);

    if(total_violation > 10. * SLEQP_ABS(func_val))
    {
      sleqp_log_warn("Problem is badly scaled, constraint violation %g significantly exceeds function value of %g",
                     total_violation,
                     func_val);
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

  if(solver->status == SLEQP_STATUS_INFEASIBLE)
  {
    sleqp_log_debug("Problem is infeasible, aborting");
    return SLEQP_OKAY;
  }

  solver->status = SLEQP_STATUS_RUNNING;

  bool reject_initial;

  SLEQP_CALL(sleqp_set_and_evaluate(problem,
                                    iterate,
                                    SLEQP_VALUE_REASON_INIT,
                                    &reject_initial));

  if(reject_initial)
  {
    sleqp_log_error("Function rejected initial solution");
    return SLEQP_INTERNAL_ERROR;
  }

  {
    SleqpSparseMatrix* cons_jac = sleqp_iterate_get_cons_jac(iterate);

    sleqp_log_info("Solving a problem with %d variables, %d constraints, %d Jacobian nonzeros",
                   num_variables,
                   num_constraints,
                   sleqp_sparse_matrix_get_nnz(cons_jac));
  }

  // ensure that the unscaled iterate is initialized
  if(solver->scaling_data)
  {
    SLEQP_CALL(sleqp_iterate_copy(iterate,
                                  solver->original_iterate));

    SleqpFunc* func = sleqp_problem_func(problem);
    const bool lsq = sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_LSQ;

    SLEQP_CALL(sleqp_unscale_iterate(solver->scaling_data,
                                     solver->original_iterate,
                                     lsq));
  }

  // Warnings
  SLEQP_CALL(print_warning(solver));

  SLEQP_CALL(sleqp_timer_reset(solver->elapsed_timer));

  solver->status = SLEQP_STATUS_RUNNING;

  solver->time_limit = time_limit;
  solver->abort_next = false;

  solver->iteration = 0;
  solver->elapsed_seconds = 0.;
  solver->last_step_type = SLEQP_STEPTYPE_NONE;

  const double deadpoint_bound = sleqp_params_get(solver->params,
                                                  SLEQP_PARAM_DEADPOINT_BOUND);

  SLEQP_CALL(sleqp_solver_print_header(solver));

  // main solving loop
  while(true)
  {
    if(time_limit != SLEQP_NONE &&
       solver->elapsed_seconds >= time_limit)
    {
      sleqp_log_info("Exhausted time limit, terminating");
      solver->status = SLEQP_STATUS_ABORT_TIME;
      break;
    }

    if(max_num_iterations != SLEQP_NONE &&
       solver->iteration >= max_num_iterations)
    {
      sleqp_log_info("Reached iteration limit, terminating");
      solver->status = SLEQP_STATUS_ABORT_ITER;
      break;
    }

    if(solver->abort_next)
    {
      sleqp_log_info("Abortion requested, terminating");
      solver->status = SLEQP_STATUS_ABORT_MANUAL;
      break;
    }

    SLEQP_CALL(sleqp_timer_start(solver->elapsed_timer));

    SLEQP_RETCODE solver_status = sleqp_solver_perform_iteration(solver);

    SLEQP_CALL(sleqp_timer_stop(solver->elapsed_timer));

    if(solver_status == SLEQP_ABORT_TIME)
    {
      sleqp_log_info("Exhausted time limit, terminating");
      solver->status = SLEQP_STATUS_ABORT_TIME;
      break;
    }

    SLEQP_CALL(solver_status);

    solver->elapsed_seconds = sleqp_timer_get_ttl(solver->elapsed_timer);

    if(solver->lp_trust_radius <= deadpoint_bound ||
       solver->trust_radius <= deadpoint_bound)
    {
      sleqp_log_warn("Reached dead point");
      solver->status = SLEQP_STATUS_ABORT_DEADPOINT;
      break;
    }

    if(solver->status != SLEQP_STATUS_RUNNING)
    {
      break;
    }
  }

  assert(solver->status != SLEQP_STATUS_RUNNING);

  double violation;

  SLEQP_CALL(sleqp_iterate_feasibility_residuum(solver->problem,
                                                solver->iterate,
                                                &violation));

  {
    SLEQP_POLISHING_TYPE polishing_type = sleqp_options_get_int(solver->options,
                                                                SLEQP_OPTION_INT_POLISHING_TYPE);

    SLEQP_CALL(sleqp_polishing_polish(solver->polishing,
                                      solver->iterate,
                                      polishing_type));
  }

  SLEQP_CALLBACK_HANDLER_EXECUTE(solver->callback_handlers[SLEQP_SOLVER_EVENT_FINISHED],
                                 SLEQP_FINISHED,
                                 solver,
                                 solver->original_iterate);

  SLEQP_CALL(sleqp_solver_print_stats(solver, violation));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solver_get_solution(SleqpSolver* solver,
                                        SleqpIterate** iterate)
{
  SLEQP_CALL(sleqp_solver_restore_original_iterate(solver));

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

  if(solver->quasi_newton)
  {
    SLEQP_CALL(sleqp_quasi_newton_reset(solver->quasi_newton));
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

static SLEQP_RETCODE solver_free(SleqpSolver** star)
{
  SleqpSolver* solver = *star;

  if(!solver)
  {
    return SLEQP_OKAY;
  }

  sleqp_free(&solver->dense_cache);

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->soc_step));

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

  SLEQP_CALL(sleqp_gauss_newton_solver_release(&solver->gauss_newton_solver));

  SLEQP_CALL(sleqp_polishing_release(&solver->polishing));

  SLEQP_CALL(sleqp_linesearch_release(&solver->linesearch));

  SLEQP_CALL(sleqp_merit_release(&solver->merit));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->estimation_residuals));

  SLEQP_CALL(sleqp_dual_estimation_free(&solver->estimation_data));

  SLEQP_CALL(sleqp_aug_jac_release(&solver->aug_jac));

  SLEQP_CALL(sleqp_sparse_factorization_release(&solver->factorization));

  SLEQP_CALL(sleqp_iterate_release(&solver->trial_iterate));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->initial_trial_point));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->trial_step));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->newton_hessian_step));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->newton_step));

  SLEQP_CALL(sleqp_newton_data_release(&solver->newton_data));

  SLEQP_CALL(sleqp_working_step_release(&solver->working_step));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->multipliers));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->cauchy_hessian_step));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->cauchy_step));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->cauchy_direction));

  SLEQP_CALL(sleqp_cauchy_release(&solver->cauchy_data));

  SLEQP_CALL(sleqp_lpi_release(&solver->lp_interface));

  SLEQP_CALL(sleqp_sparse_vector_free(&solver->original_violation));

  SLEQP_CALL(sleqp_iterate_release(&solver->iterate));

  SLEQP_CALL(sleqp_step_rule_release(&solver->step_rule));

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

  SLEQP_CALL(sleqp_quasi_newton_release(&solver->quasi_newton));

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
